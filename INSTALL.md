# INSTALL

Two ways to run this toolkit. The first takes under a minute and needs nothing but Python; the second
adds a container so that the environment is not a variable.

## 1. Audit layer only (no dependencies, no network)

Requires Python ≥ 3.9. The audit layer (`np_contracts`, `np_scm0061`, `np_report`, `npcore`) imports
only the standard library — there is nothing to install.

```bash
git clone <REPOSITORY_URL>            # or unzip the archived release
cd network-pharmacology-pro
python -m pytest tests -q             # expect: 121 passed, 1 skipped
bash demo/run_demo.sh                 # expect: DEMO OK
```

`run_demo.sh` works on a temporary copy of `demo/case`, so the bundled case is never modified and the
archive keeps passing its per-file hash check.

Expected demo output (three checkable values):

```
[contracts] PASS：26 个契约，存在 21，硬失败 0
  条目 1185 | 等级分布 {'L1': 281, 'L2': 904}
SCM0061 评分：26/34 = 76.5%  等级 B
== DEMO OK ==
```

## 2. Container (recommended for third-party verification)

```bash
docker build -f Dockerfile.slim --build-arg BASE_IMAGE=python:3.13-slim -t np-pro-slim .
docker run --rm np-pro-slim           # runs run_demo.sh, prints DEMO OK
```

The build fails if the audit layer does not pass, so a successful build is itself a check. No pip
package is installed; image size is dominated by the base image (~180 MB).

## 3. Full pipeline (network access required)

The twelve-stage pipeline queries public services (TCMSP, STRING, Open Targets, Enrichr) and needs
third-party libraries:

```bash
python -m pip install -r requirements.txt
python scripts/np_pipeline.py --inputs <inputs.json> -o <outdir>
```

Optional accelerators, detected at runtime and not required: `openpyxl`, `networkx`, `matplotlib`,
`rdkit`, AutoDock Vina (docking stage), GROMACS or OpenMM (MD stage). Missing tools cause the
corresponding stage to be skipped with a recorded reason rather than to fail silently.

## 4. Versions actually used for the reported results

| Component | Version |
|---|---|
| Python | 3.13.14 |
| OS | Windows 11 (10.0.26200) |
| Container base | `python:3.13-slim` |
| STRING (PPI) | v12.5 |
| pytest | 8.x, 121 passed / 1 skipped |
| ruff (error-level rules, `ruff.toml`) | 0 findings |

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `docker build` cannot pull the base image | Docker Hub unreachable from some networks | pass `--build-arg BASE_IMAGE=<mirror>/library/python:3.13-slim` |
| demo reports a contract hard failure | the case directory was modified by a previous manual run | use a fresh copy; `run_demo.sh` does this automatically |
| a pipeline stage is skipped | optional tool absent | see the stage's `run_manifest.json` entry for the recorded reason |
