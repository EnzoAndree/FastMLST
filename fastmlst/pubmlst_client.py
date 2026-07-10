"""Resilient authenticated HTTP client for the PubMLST REST API.

PubMLST uses short-lived OAuth1 session tokens for API requests.  This module
keeps acquisition of those tokens separate from request handling: callers
provide the long-lived access credentials and a callback that exchanges them
for a fresh session token.
"""

from __future__ import absolute_import

import json
import threading
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from numbers import Real
from urllib.parse import urlparse

import requests
from oauthlib.oauth1 import SIGNATURE_TYPE_QUERY
from requests_oauthlib import OAuth1Session


DEFAULT_TIMEOUT = 120
DEFAULT_MAX_RETRIES = 3
DEFAULT_MAX_RETRY_DELAY = 120
DEFAULT_USER_AGENT = "FastMLST"
DEFAULT_ALLOWED_HOSTS = ("rest.pubmlst.org",)


class PubMLSTClientError(RuntimeError):
    """Base class for PubMLST client failures."""


class PubMLSTAuthenticationError(PubMLSTClientError):
    """OAuth credentials or a PubMLST session token could not be used."""


class PubMLSTRequestError(PubMLSTClientError):
    """The request could not be completed at the transport layer."""


class PubMLSTHTTPError(PubMLSTClientError):
    """PubMLST returned an unsuccessful HTTP response."""

    def __init__(self, message, status_code, url):
        super().__init__(message)
        self.status_code = status_code
        self.url = url


class PubMLSTResponseError(PubMLSTClientError):
    """PubMLST returned a response that could not be decoded."""


def exponential_backoff(attempt):
    """Return the default delay for a one-based transient retry number."""
    return 2 ** (attempt - 1)


class PubMLSTClient:
    """OAuth1 PubMLST client with session renewal and bounded retries.

    ``session_token_provider`` is called with an OAuth1 session signed with the
    supplied long-lived access token.  It must return either
    ``(session_token, session_secret)`` or a mapping containing the OAuth keys
    ``oauth_token`` and ``oauth_token_secret``.

    ``session_factory``, ``backoff`` and ``sleep`` are injectable so request
    behaviour can be tested without network access or real delays.
    """

    def __init__(
        self,
        client_id,
        client_secret,
        access_token,
        access_secret,
        session_token_provider,
        timeout=DEFAULT_TIMEOUT,
        max_retries=DEFAULT_MAX_RETRIES,
        backoff=None,
        sleep=None,
        session_factory=None,
        max_retry_delay=DEFAULT_MAX_RETRY_DELAY,
        user_agent=DEFAULT_USER_AGENT,
        allowed_hosts=DEFAULT_ALLOWED_HOSTS,
    ):
        self.client_id = self._required_text(client_id, "client_id")
        self.client_secret = self._required_text(client_secret, "client_secret")
        self.access_token = self._required_text(access_token, "access_token")
        self.access_secret = self._required_text(access_secret, "access_secret")

        if not callable(session_token_provider):
            raise ValueError("session_token_provider must be callable.")
        if isinstance(max_retries, bool) or not isinstance(max_retries, int):
            raise ValueError("max_retries must be a non-negative integer.")
        if max_retries < 0:
            raise ValueError("max_retries must be a non-negative integer.")
        self._validate_timeout(timeout)
        if max_retry_delay is not None:
            if (
                isinstance(max_retry_delay, bool)
                or not isinstance(max_retry_delay, Real)
                or max_retry_delay < 0
            ):
                raise ValueError("max_retry_delay must be non-negative or None.")

        self.timeout = timeout
        self.max_retries = max_retries
        self.max_retry_delay = max_retry_delay
        self.user_agent = self._required_text(user_agent, "user_agent")
        if not allowed_hosts:
            raise ValueError("allowed_hosts must contain at least one hostname.")
        self.allowed_hosts = frozenset(
            self._required_text(host, "allowed host").lower()
            for host in allowed_hosts
        )
        self._session_token_provider = session_token_provider
        self._session_factory = session_factory or OAuth1Session
        self._backoff = backoff or exponential_backoff
        self._sleep = sleep or time.sleep
        if not callable(self._session_factory):
            raise ValueError("session_factory must be callable.")
        if not callable(self._backoff):
            raise ValueError("backoff must be callable.")
        if not callable(self._sleep):
            raise ValueError("sleep must be callable.")

        self._session = None
        self._session_generation = 0
        self._session_lock = threading.RLock()

    @staticmethod
    def _required_text(value, name):
        if not isinstance(value, str) or not value.strip():
            raise ValueError("{} is required.".format(name))
        return value

    @staticmethod
    def _validate_timeout(timeout):
        if isinstance(timeout, tuple):
            if len(timeout) != 2:
                raise ValueError("timeout tuple must contain connect and read values.")
            values = timeout
        else:
            values = (timeout,)
        for value in values:
            if isinstance(value, bool) or not isinstance(value, Real) or value <= 0:
                raise ValueError("timeout values must be positive numbers.")

    def _new_oauth_session(self, owner_token, owner_secret):
        try:
            session = self._session_factory(
                client_key=self.client_id,
                client_secret=self.client_secret,
                resource_owner_key=owner_token,
                resource_owner_secret=owner_secret,
                signature_type=SIGNATURE_TYPE_QUERY,
            )
            session.headers.update({"User-Agent": self.user_agent})
        except Exception as exc:
            raise PubMLSTAuthenticationError(
                "Unable to create the PubMLST OAuth1 session."
            ) from exc
        return session

    @staticmethod
    def _parse_session_tokens(value):
        if isinstance(value, dict):
            token = value.get("oauth_token")
            secret = value.get("oauth_token_secret")
        elif isinstance(value, (tuple, list)) and len(value) == 2:
            token, secret = value
        else:
            token, secret = None, None
        if not isinstance(token, str) or not token.strip():
            raise PubMLSTAuthenticationError(
                "The PubMLST session token provider returned no oauth_token."
            )
        if not isinstance(secret, str) or not secret.strip():
            raise PubMLSTAuthenticationError(
                "The PubMLST session token provider returned no oauth_token_secret."
            )
        return token, secret

    @staticmethod
    def _close_session(session):
        close = getattr(session, "close", None)
        if callable(close):
            close()

    @staticmethod
    def _close_response(response):
        close = getattr(response, "close", None)
        if callable(close):
            close()

    def _renew_session(self, expected_generation=None):
        """Acquire and install a fresh resource session.

        When multiple requests observe the same expired session, only the first
        one performs an exchange.  The others reuse the newly installed session.
        """
        with self._session_lock:
            if (
                expected_generation is not None
                and expected_generation != self._session_generation
                and self._session is not None
            ):
                return

            access_session = self._new_oauth_session(
                self.access_token, self.access_secret
            )
            try:
                try:
                    value = self._session_token_provider(access_session)
                except PubMLSTClientError:
                    raise
                except Exception as exc:
                    raise PubMLSTAuthenticationError(
                        "Unable to obtain a PubMLST OAuth session token."
                    ) from exc
                session_token, session_secret = self._parse_session_tokens(value)
                new_session = self._new_oauth_session(
                    session_token, session_secret
                )
            finally:
                self._close_session(access_session)

            old_session = self._session
            self._session = new_session
            self._session_generation += 1
            if old_session is not None:
                self._close_session(old_session)

    def _current_session(self):
        with self._session_lock:
            if self._session is None:
                self._renew_session()
            return self._session, self._session_generation

    def _backoff_delay(self, retry_number):
        try:
            delay = self._backoff(retry_number)
        except Exception as exc:
            raise PubMLSTClientError(
                "The PubMLST retry backoff callback failed."
            ) from exc
        if isinstance(delay, bool) or not isinstance(delay, Real) or delay < 0:
            raise PubMLSTClientError(
                "The PubMLST retry backoff callback must return a non-negative number."
            )
        return float(delay)

    def _retry_after_delay(self, response):
        headers = getattr(response, "headers", {}) or {}
        value = headers.get("Retry-After")
        if value is None:
            return None
        value = str(value).strip()
        try:
            seconds = float(value)
            if seconds < 0:
                return None
            return seconds
        except ValueError:
            pass

        try:
            retry_at = parsedate_to_datetime(value)
            if retry_at.tzinfo is None:
                retry_at = retry_at.replace(tzinfo=timezone.utc)
            return max(0.0, (retry_at - datetime.now(timezone.utc)).total_seconds())
        except (TypeError, ValueError, OverflowError):
            return None

    def _retry_delay(self, response, retry_number):
        delay = None
        if response.status_code == 429:
            delay = self._retry_after_delay(response)
        if delay is None:
            delay = self._backoff_delay(retry_number)
        if self.max_retry_delay is not None:
            delay = min(delay, float(self.max_retry_delay))
        return delay

    @staticmethod
    def _response_excerpt(response):
        try:
            text = (response.text or "").strip().replace("\n", " ")
        except Exception:
            return ""
        if not text:
            return ""
        if len(text) > 240:
            text = text[:237] + "..."
        return ": {}".format(text)

    def _http_error(self, response, url, retry_count):
        status = response.status_code
        if status in (429,) or 500 <= status <= 599:
            message = (
                "PubMLST request failed after {} transient retries "
                "(HTTP {}) for {}{}"
            ).format(
                retry_count,
                status,
                url,
                self._response_excerpt(response),
            )
        else:
            message = "PubMLST request failed (HTTP {}) for {}{}".format(
                status, url, self._response_excerpt(response)
            )
        return PubMLSTHTTPError(message, status, url)

    def _request(self, url, accept):
        if not isinstance(url, str) or not url.strip():
            raise ValueError("url is required.")
        parsed = urlparse(url)
        if parsed.scheme.lower() != "https" or (parsed.hostname or "").lower() not in self.allowed_hosts:
            raise ValueError(
                "Refusing to send PubMLST OAuth credentials to an untrusted URL: {}"
                .format(url)
            )

        transient_retries = 0
        renewed_after_401 = False
        while True:
            session, generation = self._current_session()
            try:
                response = session.get(
                    url,
                    headers={"Accept": accept},
                    timeout=self.timeout,
                    allow_redirects=False,
                )
            except (requests.RequestException, OSError) as exc:
                if transient_retries >= self.max_retries:
                    raise PubMLSTRequestError(
                        "Unable to contact PubMLST after {} transient retries "
                        "at {} ({}).".format(
                            transient_retries, url, type(exc).__name__
                        )
                    ) from None
                transient_retries += 1
                self._sleep(self._backoff_delay(transient_retries))
                continue

            try:
                status = int(response.status_code)
            except (AttributeError, TypeError, ValueError) as exc:
                raise PubMLSTResponseError(
                    "PubMLST returned a response without a valid HTTP status for {}."
                    .format(url)
                ) from exc

            if status == 401:
                if renewed_after_401:
                    self._close_response(response)
                    raise PubMLSTAuthenticationError(
                        "PubMLST rejected the renewed OAuth session (HTTP 401) for {}."
                        .format(url)
                    )
                self._close_response(response)
                self._renew_session(expected_generation=generation)
                renewed_after_401 = True
                continue

            if status == 429 or 500 <= status <= 599:
                if transient_retries >= self.max_retries:
                    error = self._http_error(response, url, transient_retries)
                    self._close_response(response)
                    raise error
                transient_retries += 1
                delay = self._retry_delay(response, transient_retries)
                self._close_response(response)
                self._sleep(delay)
                continue

            if status >= 300:
                error = self._http_error(response, url, transient_retries)
                self._close_response(response)
                raise error
            return response

    def get_json(self, url):
        """GET ``url`` and decode an authenticated JSON response."""
        response = self._request(url, "application/json")
        try:
            try:
                return response.json()
            except (ValueError, json.JSONDecodeError) as exc:
                raise PubMLSTResponseError(
                    "PubMLST returned invalid JSON for {}{}.".format(
                        url, self._response_excerpt(response)
                    )
                ) from exc
        finally:
            self._close_response(response)

    def get_text(self, url, content_type="text/plain"):
        """GET ``url`` and return its text using ``content_type`` as Accept."""
        if not isinstance(content_type, str) or not content_type.strip():
            raise ValueError("content_type is required.")
        response = self._request(url, content_type)
        try:
            try:
                return response.text
            except Exception as exc:
                raise PubMLSTResponseError(
                    "PubMLST returned an unreadable text response for {}.".format(url)
                ) from exc
        finally:
            self._close_response(response)

    def ensure_authenticated(self):
        """Obtain an OAuth session token now, without requesting a resource."""
        self._current_session()
        return self

    def close(self):
        """Close the current HTTP session, if one has been created."""
        with self._session_lock:
            session = self._session
            self._session = None
            if session is not None:
                self._close_session(session)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
