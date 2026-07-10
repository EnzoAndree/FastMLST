# Clostridioides difficile R20291 integration fixture

This fixture exercises FastMLST end to end with real, versioned biological
data and no network access during the test.

- Genome: *Clostridioides difficile* R20291 complete genome, NCBI GenBank
  accession [`FN545816.1`](https://www.ncbi.nlm.nih.gov/nuccore/FN545816.1),
  4,191,339 bases. The fixture is stored as deterministic gzip data so the test
  also exercises compressed FASTA input.
- Uncompressed genome SHA-256:
  `4afb74a527e97a21b09f8e00945db177d1693bd61424b9e5a578d452d5c8a003`.
- Compressed fixture SHA-256:
  `51c5922aad07d8ddb38c882a1297c657b24bd6b5fbc2445ebee525348f2b07e8`.
- Scheme: PubMLST `pubmlst_cdifficile_seqdef:1` (`MLST`), downloaded on
  2026-07-10. The snapshot contains 1,314 profiles and seven loci.

Expected FastMLST result:

```text
Clostridioides_difficile_R20291_FN545816.1.fasta.gz,pubmlst_cdifficile_seqdef_1,1,adk(1),atpA(1),dxr(1),glyA(10),recA(1),sodA(3),tpi(5),mlst_clade(2)
```

The integration test copies the scheme into a temporary database, builds the
BLAST database with the locally installed NCBI BLAST+ tools, and invokes the
FastMLST CLI against the genome fixture.
