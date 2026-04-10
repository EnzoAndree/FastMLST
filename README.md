![Conda](https://img.shields.io/conda/pn/bioconda/fastmlst)![CircleCI](https://img.shields.io/circleci/build/github/EnzoAndree/FastMLST/master)![GitHub](https://img.shields.io/github/license/EnzoAndree/FastMLST)[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/fastmlst/README.html)[![Downloads](https://img.shields.io/conda/dn/bioconda/fastmlst.svg?style=flat)](http://bioconda.github.io/recipes/fastmlst/README.html)![Citations](https://img.shields.io/badge/citations-44-blue)

# FastMLST
A multi-core Multilocus Sequence Typing tool coupled with allele concatenation.
# Introduction
FastMLST is a high speed standalone script wrote in Python3, which takes assemblies in FASTA format (gzipped is also allowed) and determines its ST according to MLST schemes defined in [PubMLST](https://doi.org/10.12688/wellcomeopenres.14826.1). The main advantage over other ST determination programs is that FastMLST allows the generation of a FASTA file containing the concatenated alleles for all analyzed genomes ready to be aligned and used in phylogenetic inference.

You can read a complete guide to MLST analysis in our [Wiki](https://github.com/EnzoAndree/FastMLST/wiki/The-definitive-guide-to-MLST-analysis).

## PubMLST data model and the REST API (current release)

FastMLST no longer relies on legacy static PubMLST bundle files. It uses the official **[BIGSdb REST API](https://rest.pubmlst.org)** (`https://rest.pubmlst.org`) for catalog discovery, scheme metadata, profiles, and allele FASTA downloads.

### What ships inside the package

- A **snapshot catalog** is included under `fastmlst/bundle/` (`scheme_catalog.json` and `scheme_catalog_meta.json`). This lets you run **`--scheme-list`** without network access when that snapshot is present.
- The snapshot is a **point-in-time list** of schemes; it is not a substitute for a live database. **Refreshing the catalog** (`--scheme-list-update`) and **downloading schemes** (`--update-mlst …`) require contacting the API.

### What lives in your cache

- By default, databases and downloads go to **`~/.cache/fastmlst/pubmlst`** (or **`--db_path`**).
- On a full **`--update-mlst ALL`**, the local **`scheme_catalog.json`** (and meta file) are **preserved** when the rest of the tree is reset, so listing and updates do not unnecessarily re-crawl the full API catalog every time.
- Each scheme is stored under a **stable folder name**: **`{database}_{scheme_id}`**, e.g. `pubmlst_cdifficile_seqdef_1`. Use that name with **`--scheme`**.

### When network / API access is required

| Action | Needs live API? |
|--------|-----------------|
| **`--scheme-list`** | Only if there is **no** usable catalog in cache **and** none in the installed package bundle (first run or minimal install). |
| **`--scheme-list-update`** | **Yes** — always refreshes from the API. |
| **`--update-mlst VALUE`** | **Yes** — pass **`ALL`** (downloads every scheme in the catalog; very slow) or explicit comma-separated selectors (`database:scheme_id` / stable codenames). |
| Typing genomes (after DB is built) | **No** — uses local BLAST DB and scheme files. |

### Anonymous access vs OAuth

FastMLST can query the PubMLST API without logging in, but that **anonymous mode should not be used for normal operation**. It is incomplete and may return outdated or missing data. PubMLST may limit or deny some resources to unauthenticated requests, so profile tables or allele FASTA downloads for some databases can fail or come back incomplete.

In practice:

- Use **OAuth** for **`--scheme-list-update`** and **`--update-mlst`**.
- If you need current PubMLST data, including **post-2024 allele data**, **OAuth is required**.
- The bundled or cached scheme catalog can still be viewed locally without network access, but any live refresh or download should be done with OAuth.

OAuth setup is a **one-time step**:

```bash
fastmlst --pubmlst-connect \
  --pubmlst-client-id YOUR_CLIENT_ID \
  --pubmlst-client-secret YOUR_CLIENT_SECRET
```

You can also provide the same values through **`FASTMLST_PUBMLST_CLIENT_ID`** and **`FASTMLST_PUBMLST_CLIENT_SECRET`**.

FastMLST will open the PubMLST authorization flow in your browser. After you approve access, the returned tokens are stored locally and reused in later runs, so you do not need to repeat the login every time. See **`fastmlst --help`** for the full set of OAuth-related flags.

For a step-by-step walkthrough with screenshots, see [PubMLST-API-setup-tutorial-for-FastMLST](https://github.com/EnzoAndree/FastMLST/wiki/PubMLST-API-setup-tutorial-for-FastMLST).

### Other data notes

- **`pubmlst_rmlst_seqdef`** schemes are **not** merged into the shared **`mlst.fasta`** BLAST database (rMLST loci are kept separate on purpose). They remain on disk if you install them; standard MLST BLAST typing uses the concatenated non-rMLST panel.
- Use **`--installed-scheme-stats`** for a one-line summary per installed scheme (description, inferred type, ST and allele counts).

# Installation
You can install FastMLST using either Conda or pip. If you want the most up-to-date version, you can install it directly from GitHub using pip.

### Using Conda
```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install fastmlst
```

### Using pip for the latest version
To install the latest version directly from GitHub, use the following command:
```bash
pip install git+https://github.com/EnzoAndree/FastMLST.git
```
## Dependencies
It is expected that all dependencies will be resolved when using conda for installation.
* Python > 3
* Biopython
* tqdm
* pandas
* **requests** and **requests-oauthlib** (PubMLST REST API and OAuth)
* NCBI BLAST+
* **Network access** whenever the REST API must be contacted (see table above)
# Quick Start
Examples below use short scheme labels such as `cdifficile` for readability. After **`--update-mlst …`**, the auto-detected or explicit **`--scheme`** value is typically the **stable codename** (e.g. `pubmlst_cdifficile_seqdef_1`). Run **`fastmlst --scheme-list`** to see the exact name for your install.

```
$ fastmlst cdiff_refferences/RT078_CDM120.fasta
RT078_CDM120.fasta,cdifficile,11,adk(5),atpA(8),dxr(5),glyA(11),recA(9),sodA(11),tpi(8),mlst_clade(5.0)

$ fastmlst cdiff_refferences/RT078_CDM120.fasta.gz
RT078_CDM120.fasta.gz,cdifficile,11,adk(5),atpA(8),dxr(5),glyA(11),recA(9),sodA(11),tpi(8),mlst_clade(5.0)

$ fastmlst cdiff_refferences/*.fasta
RT001_BI9.fasta,cdifficile,3,adk(1),atpA(1),dxr(2),glyA(1),recA(1),sodA(1),tpi(1),mlst_clade(1.0)
RT001_Liv24.fasta,cdifficile,3,adk(1),atpA(1),dxr(2),glyA(1),recA(1),sodA(1),tpi(1),mlst_clade(1.0)
RT002_TL178.fasta,cdifficile,8,adk(1),atpA(1),dxr(2),glyA(6),recA(1),sodA(5),tpi(1),mlst_clade(1.0)
RT012_CD630_chr_V12.fasta,cdifficile,54,adk(1),atpA(4),dxr(7),glyA(1),recA(1),sodA(3),tpi(3),mlst_clade(1.0)
RT014_TL176_v3.fasta,cdifficile,13,adk(1),atpA(1),dxr(6),glyA(1),recA(5),sodA(3),tpi(1),mlst_clade(1.0)
RT015_TL174.fasta,cdifficile,44,adk(2),atpA(5),dxr(2),glyA(1),recA(1),sodA(3),tpi(1),mlst_clade(1.0)
RT017_CF5.fasta,cdifficile,86,adk(3),atpA(7),dxr(3),glyA(8),recA(6),sodA(19),tpi(11),mlst_clade(4.0)
RT017_M68.fasta,cdifficile,37,adk(3),atpA(7),dxr(3),glyA(8),recA(6),sodA(9),tpi(11),mlst_clade(4.0)
RT023_CD305.fasta,cdifficile,791,adk(65),atpA(1),dxr(4),glyA(7),recA(2),sodA(8),tpi(7),mlst_clade(nan)
RT027_CD196.fasta,cdifficile,1,adk(1),atpA(1),dxr(1),glyA(10),recA(1),sodA(3),tpi(5),mlst_clade(2.0)
RT027_R20291_July2013.fasta,cdifficile,1,adk(1),atpA(1),dxr(1),glyA(10),recA(1),sodA(3),tpi(5),mlst_clade(2.0)
RT078_CDM120.fasta,cdifficile,11,adk(5),atpA(8),dxr(5),glyA(11),recA(9),sodA(11),tpi(8),mlst_clade(5.0)
RT106_Liv22.fasta,cdifficile,42,adk(1),atpA(1),dxr(2),glyA(1),recA(1),sodA(7),tpi(1),mlst_clade(1.0)
```
# Usage
FastMLST uses as input a assembly in FASTA format. Optionally it can be compressed with gzip or bzip2.
```
$ fastmlst cdiff_refferences/RT078_CDM120.fasta
RT078_CDM120.fasta,cdifficile,11,adk(5),atpA(8),dxr(5),glyA(11),recA(9),sodA(11),tpi(8),mlst_clade(5.0)
```
The output is a comma separated file (csv) by default, but it can be modified using the `-s` option.
```
$ fastmlst -s '\t' cdiff_refferences/RT078_CDM120.fasta
RT078_CDM120.fasta      cdifficile      11      adk(5)  atpA(8) dxr(5)  glyA(11)        recA(9) sodA(11)        tpi(8)  mlst_clade(5.0)
```
There are two options for saving the result in a text file:
```
$ fastmlst -to mlst.csv cdiff_refferences/RT078_CDM120.fasta
$ fastmlst cdiff_refferences/RT078_CDM120.fasta > mlst.csv
```
Both options generate the `mlst.csv` file containing the FastMLST result.

FastMLST is able to generate a file in FASTA format with the alleles concatenated in the same way as they are in PubMLST. If any genome is not found in this result, it means that (1) Allele contain Ns, (2) alleles missing or (3) contamination (multiple alleles for one genome). Optionally the name could be modified with `-fo` option:
```
$ fastmlst cdiff_refferences/RT078_CDM120.fasta
```
FastMLST will  try to use all available cores. It can be modified with `-t` option:
```
$ fastmlst -t 2 cdiff_refferences/RT078_CDM120.fasta 
```
You also can specify to FastMLST the scheme name using the `--scheme` option, this is particularly useful when there is more than one scheme per species. If you use this option, it will generate a table with a new format (available since version 0.0.10) which is easier to use in other programs like [phyloviz](http://www.phyloviz.net/).

```
$ fastmlst --scheme pubmlst_cdifficile_seqdef_1 cdiff_refferences/RT078_CDM120.fasta
Genome,Scheme,ST,adk,atpA,dxr,glyA,recA,sodA,tpi,mlst_clade
RT078_CDM120.fasta,cdifficile,11,5,8,5,11,9,11,8,5.0
```

If you want the old format just add the option `--legacy`:

```
$ fastmlst --legacy --scheme pubmlst_cdifficile_seqdef_1 cdiff_refferences/RT078_CDM120.fasta
RT078_CDM120.fasta,cdifficile,11,adk(5),atpA(8),dxr(5),glyA(11),recA(9),sodA(11),tpi(8),mlst_clade(5.0)
```

PubMLST schemes are listed with **`--scheme-list`** (uses your local `scheme_catalog.json` cache, the catalog shipped with the package, or a one-time API fetch if nothing is available yet). Each line includes the stable codename, `database:scheme_id`, species, and description. To **rebuild the catalog from the live API**, use **`--scheme-list-update`** (OAuth recommended; see `--pubmlst-connect`).

**Hint: use the scheme folder name (`codename`, e.g. `pubmlst_cdifficile_seqdef_1`) with `--scheme`.**

```
$ fastmlst --scheme-list
Total remote schemes: 235 (bundled with package, snapshot 2026-04-09T16:44:59+00:00)

(1) pubmlst_achromobacter_seqdef_1 | pubmlst_achromobacter_seqdef:1 | achromobacter | MLST
(2) pubmlst_abaumannii_seqdef_1 | pubmlst_abaumannii_seqdef:1 | abaumannii | MLST (Oxford)
(3) pubmlst_abaumannii_seqdef_2 | pubmlst_abaumannii_seqdef:2 | abaumannii | MLST (Pasteur)
(4) pubmlst_abaumannii_seqdef_3 | pubmlst_abaumannii_seqdef:3 | abaumannii | cgMLST v1
(5) pubmlst_actinobacillus_seqdef_1 | pubmlst_actinobacillus_seqdef:1 | actinobacillus | MLST
(6) pubmlst_aeromonas_seqdef_1 | pubmlst_aeromonas_seqdef:1 | aeromonas | MLST
(7) pubmlst_aactinomycetemcomitans_seqdef_1 | pubmlst_aactinomycetemcomitans_seqdef:1 | aactinomycetemcomitans | MLST
(8) pubmlst_aphagocytophilum_seqdef_1 | pubmlst_aphagocytophilum_seqdef:1 | aphagocytophilum | MLST
(9) pubmlst_aphagocytophilum_seqdef_2 | pubmlst_aphagocytophilum_seqdef:2 | aphagocytophilum | ankA
(10) pubmlst_aphagocytophilum_seqdef_3 | pubmlst_aphagocytophilum_seqdef:3 | aphagocytophilum | groEL
…
```

A new option in version v0.0.14 is the possibility to obtain the alleles divided into individual FASTA files (one for each allele in the scheme), ready to be used in other programs such as MLSTest.

```
$ fastmlst --scheme pubmlst_cdifficile_seqdef_1 cdiff_refferences/*.fasta --splited-output splited_mlst
Genome,Scheme,ST,adk,atpA,dxr,glyA,recA,sodA,tpi,mlst_clade
RT001_BI9.fasta,cdifficile,3,1,1,2,1,1,1,1,1.0
RT001_Liv24.fasta,cdifficile,3,1,1,2,1,1,1,1,1.0
RT002_TL178.fasta,cdifficile,8,1,1,2,6,1,5,1,1.0
RT012_CD630_chr_V12.fasta,cdifficile,54,1,4,7,1,1,3,3,1.0
RT014_TL176_v3.fasta,cdifficile,13,1,1,6,1,5,3,1,1.0
RT015_TL174.fasta,cdifficile,44,2,5,2,1,1,3,1,1.0
RT017_CF5.fasta,cdifficile,86,3,7,3,8,6,19,11,4.0
RT017_M68.fasta,cdifficile,37,3,7,3,8,6,9,11,4.0
RT023_CD305.fasta,cdifficile,791,65,1,4,7,2,8,7,
RT027_CD196.fasta,cdifficile,1,1,1,1,10,1,3,5,2.0
RT027_R20291_July2013.fasta,cdifficile,1,1,1,1,10,1,3,5,2.0
RT078_CDM120.fasta,cdifficile,11,5,8,5,11,9,11,8,5.0
RT106_Liv22.fasta,cdifficile,42,1,1,2,1,1,7,1,1.0
```

```
$ ls splited_mlst/
adk.fasta  atpA.fasta  dxr.fasta  glyA.fasta  recA.fasta  sodA.fasta  tpi.fasta
$ cat splited_mlst/adk.fasta
>RT001_BI9.fasta adk
CATATATCAACAGGAGATATATTCAGAAAGAATATAAAAGAGGGAACAGAACTTGGAAAA
AAAGCTAAAGAATACATGGACCAAGGTTTATTAGTACCAGATGAGTTAACTGTAGGTTTA
GTTACTGATAGAATATCTCAAGAAGATTGTAAAAATGGATTTATGTTAGATGGATTTCCA
AGAAATGTAGCACAAGGAGAACATTTAGATATCTTCTTAAAAAATGCTGGTATATCACTA
GATAAAGTTGTCAATATTGAAGTTGATAAGAGTATATTAGTGTCTAGAGCAGTTGGTAGA
AGAATATGTAAGTCTTGTGGAGCTACTTACCATGTTGAGTTTAATCCTCCTAAAGTAGAA
GGTGTATGTGATGTATGCCAAGGAGAATTATATCAAAGAGCTGATGATAATGAAGAAACT
GTATCTAAGAGAATACAAGTTTATCTAGATGAAACTAAGCCTTTAGTAGATTATTATAGC
AAACAAGGTATAATAGCAGAT
...
>RT106_Liv22.fasta adk
CATATATCAACAGGAGATATATTCAGAAAGAATATAAAAGAGGGAACAGAACTTGGAAAA
AAAGCTAAAGAATACATGGACCAAGGTTTATTAGTACCAGATGAGTTAACTGTAGGTTTA
GTTACTGATAGAATATCTCAAGAAGATTGTAAAAATGGATTTATGTTAGATGGATTTCCA
AGAAATGTAGCACAAGGAGAACATTTAGATATCTTCTTAAAAAATGCTGGTATATCACTA
GATAAAGTTGTCAATATTGAAGTTGATAAGAGTATATTAGTGTCTAGAGCAGTTGGTAGA
AGAATATGTAAGTCTTGTGGAGCTACTTACCATGTTGAGTTTAATCCTCCTAAAGTAGAA
GGTGTATGTGATGTATGCCAAGGAGAATTATATCAAAGAGCTGATGATAATGAAGAAACT
GTATCTAAGAGAATACAAGTTTATCTAGATGAAACTAAGCCTTTAGTAGATTATTATAGC
AAACAAGGTATAATAGCAGAT
```

## Custom MLST Database Location

FastMLST now supports configuring a custom location for the PubMLST database. By default, the tool uses a cache directory at `~/.cache/fastmlst/pubmlst`. However, if you prefer to store the database in an alternate location (for example, on a high-performance drive or in a centralized directory), you can override this default path using the `--db_path` command-line argument.

### How It Works

When the `--db_path` option is provided, FastMLST calls a helper function (`set_pathdb`) that:
- **Overrides the default database path:** The internal global `pathdb` variable is updated to use your specified path.
- **Ensures the custom directory exists:** The directory is automatically created if it does not exist.
- **Uses the custom path for all subsequent operations:** All processes (such as fetching, updating, or reading database files) use the new path.

### Usage Example

To run FastMLST with a custom MLST database directory, simply use the `--db_path` option:

```bash
$ fastmlst --db_path /path/to/your/custom/db [other-options] genomes...
```

For instance, if you want the MLST database to reside in `/data/fastmlst_db`, run:

```bash
$ fastmlst --db_path /data/fastmlst_db cdiff_refferences/RT078_CDM120.fasta
```

### When to Use This Feature

- **Optimizing I/O Performance:** Place the database on a disk with faster read/write speeds.
- **Managing Disk Usage:** Store the database on a separate partition or drive with more available space.
- **Custom Deployment Setups:** Particularly useful in multi-user or cluster environments where centralized data management is preferred.

**Note:** Ensure that the directory you specify has proper write permissions. FastMLST will automatically create the directory (and any necessary parent directories) if they do not already exist.

## Output symbology

Symbol | Meaning | Length | Identity
---   | --- | --- | ---
`n`   | Exact intact allele                   | 100%            | 100%
`~n`  | Novel full length allele similar to n | 100%            | &ge; `-pid`
`n?`  | Partial match to known allele        | &ge; `-cov` | &ge; `-pid`
`-`   | Allele missing (or allele containing Ns) | &lt; `-cov` | &lt; `-pid`
`n,m` | Multiple alleles                     | &nbsp;          | &nbsp;
## Scoring system
FastMLST uses a scoring system to determine the scheme to be employed similar to that proposed by [Tseemann](https://github.com/tseemann/mlst). The score for a scheme with N alleles is as follows:

* +100/N points for an exact allele match _e.g._ `1`
* +70/N points for a novel allele match _e.g._ `~1`
* +20/N points for a partial allele match _e.g._ `1?`
* 0 points for a missing allele _e.g._ `-`
# Updating the Schemes
You should **always, always, always keep the PubMLST database updated**. **`--update-mlst` requires a value** and uses the live REST API to download profiles and allele FASTA (OAuth recommended).

**Typical — only the schemes you need:**

```
$ fastmlst --update-mlst "pubmlst_cdifficile_seqdef:1,pubmlst_neisseria_seqdef:1"
# or stable codenames:
$ fastmlst --update-mlst "pubmlst_cdifficile_seqdef_1"
```

Discover IDs with **`fastmlst --scheme-list`**.

**Full install — every scheme in the catalog (slow, explicit):**

```
$ fastmlst --update-mlst ALL
```

Plain **`fastmlst --update-mlst`** (without a value) is **rejected** by the CLI.

If BLAST database files are missing and you run typing without updating first, FastMLST exits with an error pointing you to **`--update-mlst …`**.

The scheme catalog is taken from **local cache** or the **bundled snapshot** when available, so listing and incremental work do not unnecessarily re-crawl the full API catalog every time.

Refresh only the **scheme list** from the API (then print it):

```
$ fastmlst --scheme-list-update
```

# Complete usage Options
```
usage: fastmlst.py [-h] [-t THREADS] [-v {0,1,2}] [-s SEPARATOR] [-sch SCHEME] [--scheme-list] [--scheme-list-update] [--installed-scheme-stats] [-fo FASTAOUTPUT]
                   [-to TABLEOUTPUT] [-cov COVERAGE] [-pid IDENTITY] [--update-mlst [UPDATE_MLST]] [--redownload-all-schemes] [-sp SPLITED_OUTPUT]
                   [--fasta2line] [--longheader] [--legacy] [-n NOVEL] [-V] [--db_path DB_PATH] [--pubmlst-client-id PUBMLST_CLIENT_ID]
                   [--pubmlst-client-secret PUBMLST_CLIENT_SECRET] [--pubmlst-connect]
                   [genomes ...]

⚡️🧬 FastMLST: A multi-core tool for multilocus sequence typing of draft genome assemblies

positional arguments:
  genomes

options:
  -h, --help            show this help message and exit
  -t, --threads THREADS
                        Number of threads to use (default 14)
  -v, --verbose {0,1,2}
                        Verbose output level choices: [0, 1, 2]
  -s, --separator SEPARATOR
                        Choose a character to use as a separator (default ",")
  -sch, --scheme SCHEME
                        Set a scheme target (I am not dumb, let me choose a scheme by myself!)
  --scheme-list         List PubMLST schemes (cache, bundled catalog, or one-time API fetch if missing)
  --scheme-list-update  Refresh scheme catalog from the PubMLST API, then print the list
  --installed-scheme-stats
                        List each installed scheme on one line: PubMLST description, type, ST count, alleles, per-locus counts
  -fo, --fastaoutput FASTAOUTPUT
                        File name of the concatenated alleles output (default "")
  -to, --tableoutput TABLEOUTPUT
                        File name of the MLST table output (default STDOUT)
  -cov, --coverage COVERAGE
                        DNA %Cov to report high quality partial allele [?] (default 99%)
  -pid, --identity IDENTITY
                        DNA %Identity of full allelle to consider 'similar' [~] (default 95%)
  --update-mlst [UPDATE_MLST]
                        Update PubMLST from the API. Pass ALL for the full catalog (very slow), or a comma-separated list of database:scheme_id / stable codenames.
  --redownload-all-schemes
                        Ignore local version metadata and re-download every scheme (default: skip schemes that match remote API metadata)
  -sp, --splited-output SPLITED_OUTPUT
                        Directory output for splited alleles (default "")
  --fasta2line          The fasta files will be in fasta2line format
  --longheader          If --longheader is invoked, the header of FASTA file contain a long description
  --legacy              If --legacy is invoked, the csv reported contain the gene name and the allele id in the row
                        [adk(1),atpA(4),dxr(7),glyA(1),recA(1),sodA(3),tpi(3)]. This option is only available when the --scheme is defined
  -n, --novel NOVEL     File name of the novel alleles
  -V, --version         Show program's version number and exit
  --db_path DB_PATH     Custom directory for MLST database (default: ~/.cache/fastmlst/pubmlst)
  --pubmlst-client-id PUBMLST_CLIENT_ID
                        PubMLST OAuth client ID (or FASTMLST_PUBMLST_CLIENT_ID)
  --pubmlst-client-secret PUBMLST_CLIENT_SECRET
                        PubMLST OAuth client secret (or FASTMLST_PUBMLST_CLIENT_SECRET)
  --pubmlst-connect     Run one-time PubMLST OAuth setup and save credentials/tokens
```

# Citation

Guerrero-Araya E, Muñoz M, Rodríguez C, Paredes-Sabja D. FastMLST: A Multi-core Tool for Multilocus Sequence Typing of Draft Genome Assemblies. Bioinform Biol Insights. 2021 Nov 27;15:11779322211059238. doi: [10.1177/11779322211059238](https://doi.org/10.1177/11779322211059238). PMID: 34866905; PMCID: [PMC8637782](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc8637782/).
