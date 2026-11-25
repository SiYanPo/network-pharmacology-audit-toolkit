# Provenance of this repository

This repository is a **byte-exact mirror** of the archived release
`network-pharmacology-audit-toolkit` v3.0.19, deposited on Zenodo.

| | |
|---|---|
| Archive | `network-pharmacology-audit-toolkit-v3.0.19.zip` |
| Size | 3.15 MiB zipped, 8.13 MiB unpacked |
| Files | 636 |
| ZIP SHA-256 | `a3b7282af9f6924a12a188886a4e9eb324b2beba853a0d22512550cdc6630aa6` |
| Per-file manifest | `ARCHIVE_MANIFEST.tsv` (636 rows: sha256, bytes, path) |
| Built | 2026-11-25 |

## Why byte-exact, and what it costs

Almost every paper's GitHub repository and its Zenodo snapshot differ in some
way: a `.gitignore` added later, line endings normalised by a contributor's
editor, a file renamed after the deposit. That is convenient and it quietly
breaks the only claim an archive is supposed to support — that a reader can
confirm the code they are running is the code that was reviewed.

Here we accept a small cost to keep the claim checkable. `.gitattributes`
sets `* -text`, which disables git's end-of-line conversion. Without it, the
460 CRLF files in the tree would be rewritten to LF on a Windows checkout and
every hash in `ARCHIVE_MANIFEST.tsv` would fail.

The cost is that line-ending changes show up as whole-file diffs. We consider
that acceptable: this repository is not where development happens, it is where
the reviewed artefact lives.

## Verify it yourself

```bash
git clone <this repository> np-toolkit && cd np-toolkit

# 1. every file against the manifest recorded at deposit time
python verify_archive.py ARCHIVE_MANIFEST.tsv .
#    expect: checked 636 files: 0 mismatched, 0 missing

# 2. the tests still pass
python -m pytest tests -q
#    expect: 121 passed, 1 skipped

# 3. the demo still reaches the reported scorecard
bash demo/run_demo.sh
#    expect: contracts PASS, SCM0061 26/34 (grade B), DEMO OK
```

If step 1 fails on a fresh clone, something changed the bytes in transit and
that is worth knowing before trusting anything downstream.

## What is not here

Third-party databases (TCMSP, STRING, Open Targets, DrugBank, Europe PMC) are
not redistributed. Each pipeline stage writes an `environment.json` recording
the database version and retrieval parameters actually used, so the inputs can
be re-fetched rather than guessed.

Build junk (`__pycache__`, `.pytest_cache`, `logs/`) is excluded from the
archive. `_ruff_cache/` is retained because it was present when the archive
hash was fixed; removing it now would silently change every downstream
checksum.

## Two files that are not in the archive

`.gitattributes` (in this repository only) and this file are repository
metadata, added after the deposit. They are the only two paths here that do
not appear in `ARCHIVE_MANIFEST.tsv`. `verify_archive.py` checks the manifest
against the tree and therefore reports no mismatch; a tool that additionally
flags unlisted files will report exactly these two, which is expected.
