# NG-STAR v2 source snapshot

Offline test snapshot downloaded from the official NG-STAR v2 website on
2026-07-10.

- Loci: `penA`, `mtrR`, `porB`, `ponA`, `gyrA`, `parC`, and `23S`.
- Alleles: `https://ngstar.canada.ca/alleles/download?lang=en&loci_name=LOCUS`
- Profiles: `https://ngstar.canada.ca/sequence_types/download?lang=en`
- Profile workbook: 8,125 profiles with the columns `Sequence Type`, `penA`,
  `mtrR`, `porB`, `ponA`, `gyrA`, `parC`, and `23S`.

The official FASTA responses are stored with deterministic gzip compression.
The test decompresses them before presenting the payloads to the downloader.
The XLSX is preserved exactly as downloaded.
