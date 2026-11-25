# Archive summary

- **Toolkit**: network-pharmacology-audit-toolkit v3.0.19
- **Files**: 636
- **Uncompressed**: 8.13 MiB
- **ZIP**: `network-pharmacology-audit-toolkit-v3.0.19.zip` — 3.15 MiB — sha256 `a3b7282af9f6924a12a188886a4e9eb324b2beba853a0d22512550cdc6630aa6`
- **TAR.GZ**: `network-pharmacology-audit-toolkit-v3.0.19.tar.gz` — 2.86 MiB — sha256 `ed273f139fe4cd18a3aa389086b0e237842965b4dc1e36e84bb97b1108ee841d`

## Excluded from the archive (build junk, not source)

| Pattern | Reason |
|---|---|
| `UsershuaweiAppDataLocalTemp*/` | directories created by an earlier path-mangling bug; contain benchmark scratch output only |
| `__pycache__/`, `*.pyc` | interpreter cache |
| `.pytest_cache/`, `.coverage` | test artefacts |
| `logs/` | local run logs (no scientific content) |
| `_rs.txt`, `.DS_Store`, `Thumbs.db` | scratch files |

## File count by top-level directory

| Path | Files |
|---|---|
| `audit` | 289 |
| `demo` | 87 |
| `scripts` | 83 |
| `references` | 40 |
| `.ruff_cache` | 23 |
| `tools` | 19 |
| `tests` | 17 |
| `data` | 12 |
| `examples` | 10 |
| `benchmark` | 8 |
| `adapters` | 7 |
| `evals` | 4 |
| `.github` | 3 |
| `assets` | 3 |
| `config` | 3 |
| `contracts` | 2 |
| `policies` | 2 |
| `ARCHIVE_README.md` | 1 |
| `CHANGELOG.md` | 1 |
| `Dockerfile` | 1 |
| `Dockerfile.slim` | 1 |
| `INSTALL.md` | 1 |
| `LICENSE` | 1 |
| `PORTABILITY.md` | 1 |
| `README.md` | 1 |
| `REPRODUCIBILITY.md` | 1 |
| `SKILL.md` | 1 |
| `app.py` | 1 |
| `audit_help.json` | 1 |
| `audit_modules.json` | 1 |
| `audit_selftest.json` | 1 |
| `coverage.json` | 1 |
| `coverage_baseline.txt` | 1 |
| `docs` | 1 |
| `requirements-dev.txt` | 1 |
| `requirements-lock.txt` | 1 |
| `requirements.txt` | 1 |
| `ruff.toml` | 1 |
| `run_benchmark.py` | 1 |
| `setup_env.sh` | 1 |
| `verify_archive.py` | 1 |

## Known extractor caveat (added after the archive was built)

Three files under `data/herb_cache/` have Chinese names (`何首乌.json`, `楮实子.json`,
`龟甲.json`). The zip stores them with the UTF-8 filename flag (bit 11) set, which Python,
7-Zip, macOS Archive Utility and the Windows Explorer shell all honour. The Info-ZIP 6.0
`unzip` build shipped with Git for Windows does **not**: it decodes entry names with the OEM
code page (CP936) and writes mojibake, after which a per-file hash check reports three false
"missing" files. Use `tar -xf` (bundled with Windows 10+, macOS and Linux) or `verify_download.py`,
which hashes inside the zip and therefore never involves an extractor's code page.

```bash
python verify_download.py network-pharmacology-audit-toolkit-v3.0.19.zip ARCHIVE_MANIFEST.tsv
# expect: checked : 636 matched, 0 mismatched, 0 missing, 0 unlisted
```

## Integrity verification (third party)

Run this **before** executing anything; the demo writes only to a temporary copy, so the
hash check below still passes after a demo run (verified).

```bash
sha256sum network-pharmacology-audit-toolkit-v3.0.19.zip        # must equal the value above
unzip -q network-pharmacology-audit-toolkit-v3.0.19.zip && cd network-pharmacology-audit-toolkit-v3.0.19
python verify_archive.py PATH/TO/ARCHIVE_MANIFEST.tsv .   # shipped as verify_archive.py
```

`verify_archive.py` (three lines of logic, reproduced here for transparency) recomputes the
SHA-256 of every file and compares it with the manifest:

```python
import csv, hashlib, sys
man, root = sys.argv[1], sys.argv[2]
bad = 0
for r in csv.DictReader(open(man, encoding='utf-8'), delimiter='\t'):
    h = hashlib.sha256(open(root + '/' + r['path'].split('/', 1)[1], 'rb').read()).hexdigest()
    if h != r['sha256']:
        bad += 1; print('MISMATCH', r['path'])
print('verified OK' if bad == 0 else f'{bad} mismatches')
```

## Functional verification (after extraction, no network)

```bash
python -m pytest tests -q          # expect 121 passed / 1 skipped
bash demo/run_demo.sh              # expect 'DEMO OK' and a 26/34 scorecard
docker build -f Dockerfile.slim --build-arg BASE_IMAGE=<mirror>/library/python:3.13-slim -t np-pro-slim .
docker run --rm np-pro-slim        # expect 'DEMO OK'
```
