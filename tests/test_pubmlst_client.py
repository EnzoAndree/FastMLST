import unittest

import requests
from oauthlib.oauth1 import SIGNATURE_TYPE_QUERY

from fastmlst.pubmlst_client import PubMLSTAuthenticationError
from fastmlst.pubmlst_client import PubMLSTClient
from fastmlst.pubmlst_client import PubMLSTHTTPError
from fastmlst.pubmlst_client import PubMLSTRequestError
from fastmlst.pubmlst_client import PubMLSTResponseError


class FakeResponse:
    def __init__(self, status_code=200, json_data=None, text="", headers=None):
        self.status_code = status_code
        self.json_data = json_data
        self.text = text
        self.headers = headers or {}
        self.closed = False

    def json(self):
        if isinstance(self.json_data, Exception):
            raise self.json_data
        return self.json_data

    def close(self):
        self.closed = True


class FakeSession:
    def __init__(self, responses=None):
        self.responses = list(responses or [])
        self.headers = {}
        self.calls = []
        self.closed = False

    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        if not self.responses:
            raise AssertionError("Unexpected request to {}".format(url))
        result = self.responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    def close(self):
        self.closed = True


class FakeSessionFactory:
    def __init__(self, resource_sessions):
        self.resource_sessions = list(resource_sessions)
        self.calls = []
        self.bootstrap_sessions = []

    def __call__(self, **kwargs):
        self.calls.append(kwargs)
        if kwargs["resource_owner_key"] == "access-token":
            session = FakeSession()
            self.bootstrap_sessions.append(session)
            return session
        if not self.resource_sessions:
            raise AssertionError("Unexpected resource session creation")
        return self.resource_sessions.pop(0)


class TokenProvider:
    def __init__(self, values):
        self.values = list(values)
        self.sessions = []

    def __call__(self, access_session):
        self.sessions.append(access_session)
        if not self.values:
            raise AssertionError("Unexpected token renewal")
        result = self.values.pop(0)
        if isinstance(result, Exception):
            raise result
        return result


def make_client(resource_sessions, token_values=None, **kwargs):
    factory = FakeSessionFactory(resource_sessions)
    provider = TokenProvider(token_values or [("session-token", "session-secret")])
    client = PubMLSTClient(
        "client-id",
        "client-secret",
        "access-token",
        "access-secret",
        provider,
        session_factory=factory,
        sleep=kwargs.pop("sleep", lambda _delay: None),
        allowed_hosts=("rest.pubmlst.org", "example.test"),
        **kwargs
    )
    return client, factory, provider


class TestPubMLSTClient(unittest.TestCase):
    def test_get_json_uses_query_signature_accept_and_timeout(self):
        resource = FakeSession([FakeResponse(json_data={"schemes": [1]})])
        client, factory, provider = make_client([resource], timeout=(4, 17))

        result = client.get_json("https://rest.pubmlst.org/db")

        self.assertEqual(result, {"schemes": [1]})
        self.assertEqual(len(provider.sessions), 1)
        self.assertEqual(len(factory.calls), 2)
        for call in factory.calls:
            self.assertEqual(call["signature_type"], SIGNATURE_TYPE_QUERY)
        self.assertEqual(
            factory.calls[0]["resource_owner_key"], "access-token"
        )
        self.assertEqual(
            factory.calls[1]["resource_owner_key"], "session-token"
        )
        self.assertEqual(resource.headers["User-Agent"], "FastMLST")
        self.assertEqual(
            resource.calls,
            [
                (
                    "https://rest.pubmlst.org/db",
                    {
                        "headers": {"Accept": "application/json"},
                        "timeout": (4, 17),
                        "allow_redirects": False,
                    },
                )
            ],
        )

    def test_get_text_sets_requested_accept_header(self):
        response = FakeResponse(text=">abc\nACGT\n")
        resource = FakeSession([response])
        client, _factory, _provider = make_client([resource])

        result = client.get_text(
            "https://rest.pubmlst.org/alleles", content_type="text/x-fasta"
        )

        self.assertEqual(result, ">abc\nACGT\n")
        self.assertEqual(
            resource.calls[0][1]["headers"], {"Accept": "text/x-fasta"}
        )
        self.assertEqual(resource.calls[0][1]["timeout"], 120)
        self.assertTrue(response.closed)

    def test_429_honours_retry_after_before_succeeding(self):
        resource = FakeSession(
            [
                FakeResponse(429, text="slow down", headers={"Retry-After": "3"}),
                FakeResponse(json_data={"ok": True}),
            ]
        )
        delays = []
        client, _factory, _provider = make_client(
            [resource], sleep=delays.append, max_retries=1
        )

        self.assertEqual(client.get_json("https://example.test/data"), {"ok": True})
        self.assertEqual(delays, [3.0])
        self.assertEqual(len(resource.calls), 2)

    def test_5xx_uses_injected_backoff(self):
        resource = FakeSession(
            [
                FakeResponse(500, text="one"),
                FakeResponse(503, text="two"),
                FakeResponse(text="done"),
            ]
        )
        backoff_calls = []
        delays = []

        def backoff(attempt):
            backoff_calls.append(attempt)
            return attempt / 4.0

        client, _factory, _provider = make_client(
            [resource],
            max_retries=2,
            backoff=backoff,
            sleep=delays.append,
        )

        self.assertEqual(client.get_text("https://example.test/data"), "done")
        self.assertEqual(backoff_calls, [1, 2])
        self.assertEqual(delays, [0.25, 0.5])

    def test_transient_retries_are_bounded(self):
        resource = FakeSession(
            [FakeResponse(503, text="unavailable") for _unused in range(3)]
        )
        client, _factory, _provider = make_client(
            [resource], max_retries=2
        )

        with self.assertRaises(PubMLSTHTTPError) as caught:
            client.get_text("https://example.test/data")

        self.assertEqual(caught.exception.status_code, 503)
        self.assertEqual(len(resource.calls), 3)
        self.assertIn("after 2 transient retries", str(caught.exception))

    def test_401_renews_session_once_and_repeats_request(self):
        expired = FakeSession([FakeResponse(401, text="expired")])
        renewed = FakeSession([FakeResponse(json_data={"ok": True})])
        client, factory, provider = make_client(
            [expired, renewed],
            token_values=[
                ("session-token", "session-secret"),
                {"oauth_token": "new-token", "oauth_token_secret": "new-secret"},
            ],
        )

        result = client.get_json("https://example.test/data")

        self.assertEqual(result, {"ok": True})
        self.assertEqual(len(provider.sessions), 2)
        self.assertEqual(len(expired.calls), 1)
        self.assertEqual(len(renewed.calls), 1)
        self.assertTrue(expired.closed)
        self.assertEqual(factory.calls[-1]["resource_owner_key"], "new-token")

    def test_second_401_raises_without_another_renewal(self):
        expired = FakeSession([FakeResponse(401)])
        still_rejected = FakeSession([FakeResponse(401)])
        client, _factory, provider = make_client(
            [expired, still_rejected],
            token_values=[
                ("session-token", "session-secret"),
                ("new-token", "new-secret"),
            ],
        )

        with self.assertRaises(PubMLSTAuthenticationError) as caught:
            client.get_json("https://example.test/data")

        self.assertEqual(len(provider.sessions), 2)
        self.assertIn("renewed OAuth session", str(caught.exception))

    def test_invalid_session_token_response_is_clear(self):
        client, _factory, _provider = make_client(
            [], token_values=[{"oauth_token": "only-token"}]
        )

        with self.assertRaises(PubMLSTAuthenticationError) as caught:
            client.get_json("https://example.test/data")

        self.assertIn("oauth_token_secret", str(caught.exception))

    def test_json_decode_failure_is_wrapped(self):
        resource = FakeSession(
            [FakeResponse(json_data=ValueError("bad json"), text="not-json")]
        )
        client, _factory, _provider = make_client([resource])

        with self.assertRaises(PubMLSTResponseError) as caught:
            client.get_json("https://example.test/data")

        self.assertIn("invalid JSON", str(caught.exception))

    def test_transport_failure_is_wrapped(self):
        resource = FakeSession(
            [
                requests.ConnectionError(
                    "failed for https://example.test/data?oauth_signature=sensitive"
                )
                for _unused in range(3)
            ]
        )
        delays = []
        client, _factory, _provider = make_client(
            [resource], max_retries=2, sleep=delays.append
        )

        with self.assertRaises(PubMLSTRequestError) as caught:
            client.get_text("https://example.test/data")

        self.assertIn("ConnectionError", str(caught.exception))
        self.assertNotIn("oauth_signature", str(caught.exception))
        self.assertNotIn("sensitive", str(caught.exception))
        self.assertEqual(len(resource.calls), 3)
        self.assertEqual(delays, [1.0, 2.0])

    def test_transport_failure_is_retried_before_succeeding(self):
        resource = FakeSession(
            [requests.exceptions.SSLError('temporary TLS EOF'), FakeResponse(text='ok')]
        )
        delays = []
        client, _factory, _provider = make_client(
            [resource], max_retries=1, sleep=delays.append
        )

        self.assertEqual(client.get_text('https://example.test/data'), 'ok')
        self.assertEqual(len(resource.calls), 2)
        self.assertEqual(delays, [1.0])

    def test_non_retryable_4xx_fails_immediately(self):
        resource = FakeSession([FakeResponse(404, text="missing")])
        client, _factory, _provider = make_client([resource], max_retries=5)

        with self.assertRaises(PubMLSTHTTPError) as caught:
            client.get_text("https://example.test/data")

        self.assertEqual(caught.exception.status_code, 404)
        self.assertEqual(len(resource.calls), 1)
        self.assertIn("missing", str(caught.exception))

    def test_close_releases_current_session(self):
        resource = FakeSession([FakeResponse(text="ok")])
        client, _factory, _provider = make_client([resource])
        client.get_text("https://example.test/data")

        client.close()

        self.assertTrue(resource.closed)

    def test_ensure_authenticated_obtains_token_without_resource_request(self):
        resource = FakeSession([FakeResponse(text="unused")])
        client, factory, provider = make_client([resource])

        result = client.ensure_authenticated()

        self.assertIs(result, client)
        self.assertEqual(len(provider.sessions), 1)
        self.assertEqual(len(factory.calls), 2)
        self.assertEqual(resource.calls, [])

    def test_configuration_is_validated(self):
        provider = TokenProvider([("token", "secret")])
        factory = FakeSessionFactory([])
        with self.assertRaises(ValueError):
            PubMLSTClient(
                "", "secret", "access", "access-secret", provider,
                session_factory=factory,
            )
        with self.assertRaises(ValueError):
            PubMLSTClient(
                "id", "secret", "access", "access-secret", provider,
                timeout=0, session_factory=factory,
            )
        with self.assertRaises(ValueError):
            PubMLSTClient(
                "id", "secret", "access", "access-secret", provider,
                max_retries=-1, session_factory=factory,
            )

    def test_rejects_untrusted_or_non_https_urls_before_signing(self):
        client, _factory, _provider = make_client([])
        with self.assertRaisesRegex(ValueError, "untrusted URL"):
            client.get_json("https://attacker.example/data")
        with self.assertRaisesRegex(ValueError, "untrusted URL"):
            client.get_json("http://rest.pubmlst.org/db")


if __name__ == "__main__":
    unittest.main()
