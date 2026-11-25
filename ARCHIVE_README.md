# network-pharmacology-audit-toolkit

An **audit-first** network-pharmacology toolkit: a 12-stage pipeline in which every intermediate
product is contract-checked, every compound–target edge carries a graded evidence level (L0–L4), and
every run is scored against the evaluation items of **SCM0061-2021** (the *Guidelines for evaluation
methods of network pharmacology*, World Federation of Chinese Medicine Societies, 2021).

This archive is the exact snapshot referenced by the manuscript submitted to *Journal of
Cheminformatics* (Software article). It is intended to be **runnable without the original author's
environment**: the audit layer depends only on the Python standard library.

## What you can verify in five minutes, offline

```bash
python -m pytest tests -q      # 121 passed, 1 skipped
bash demo/run_demo.sh          # contract check -> evidence ledger -> SCM0061 scorecard -> HTML report
```

Expected demo output (byte-stable modulo report size):

```
[contracts] PASS：26 个契约，存在 21，硬失败 0
SCM0061 评分：26/34 = 76.5%  等级 B
== DEMO OK ==
```

Container route (nothing installed but the base image):

```bash
docker build -f Dockerfile.slim --build-arg BASE_IMAGE=python:3.13-slim -t np-pro-slim .
docker run --rm np-pro-slim      # -> DEMO OK
```

## Layout

| Path | Contents |
|---|---|
| `scripts/` | Pipeline modules (`np_*.py`), including the audit layer (`np_contracts`, `np_scm0061`, `np_evidence*`) |
| `config/module_registry.yaml` | Module registry: when each module runs and its minimal command |
| `contracts/` | Product contracts (26): required columns, types, optionality |
| `demo/` | The bundled worked case (quercetin × colorectal cancer) and the offline demo script |
| `tests/` | 121 regression tests (contracts, registry, honesty guards) |
| `tools/` | Auxiliary analyses (trajectory audit, DrugBank recall audit, bibliometrics, replicate merging) |
| `references/` | Protocols: MD reuse and replicate statistics, external data sources |
| `docs/` | Code-availability statement draft |
| `benchmark/`, `data/` | Internal golden set and cached reference data |
| `audit/` | Development-time audit artefacts (cross-checks against third-party tools) |

## What this archive deliberately does **not** do

- It does not ship third-party databases (TCMSP, STRING, DrugBank): their licences forbid
  redistribution. Cached extracts under `data/` are limited to what is needed to reproduce the
  bundled case.
- It does not claim to make network-pharmacology conclusions more accurate. It makes their
  **reproducibility status diagnosable** and their **failure modes attributable**; the manuscript
  reports a null result for accuracy alongside the positive result for feasibility.

## Licence

MIT (see `LICENSE`). The copyright holder string reads "The Authors" and must be completed before
publication of the archived record.

## Integrity

`ARCHIVE_MANIFEST.tsv` (shipped next to the archive, not inside it) lists the SHA-256 of every file
together with its byte size. To verify:

```bash
sha256sum network-pharmacology-audit-toolkit-v3.0.17.zip      # compare with ARCHIVE_SUMMARY.md
unzip -q network-pharmacology-audit-toolkit-v3.0.17.zip
python network-pharmacology-audit-toolkit-v3.0.17/verify_archive.py \
       ARCHIVE_MANIFEST.tsv network-pharmacology-audit-toolkit-v3.0.17
```

Do this **before** running the demo. `demo/run_demo.sh` works on a temporary copy of the case so that
`demo/case/` is never rewritten — but any manual re-run of `scripts/np_scm0061.py --dir demo/case`
would write the scorecard into the case and change its hash (this happened during preparation and is
why the demo was changed to use a copy).
