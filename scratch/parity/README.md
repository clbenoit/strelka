# Py2/Py3 parity validation — EVS training tools

Fixed data set, fixed training chain and fixed C++ scoring runs used to
validate the Python 3 port of `scoringModelTraining/`
(change `strelka-py3-training`). The clinical acceptance criterion is parity
between the unconverted Python 2 code (reference tool chain) and the
converted Python 3 code (pinned training environment):

* **CP1** — feature CSVs byte-identical
* **CP2** — exported model JSONs: feature map identical + tree divergences
  quantified (`compare_models.py`) + per-variant P(TP) scored by each
  tool chain's own model (`score_check.py`)
* **CP3** — both model JSONs loaded in C++ (`RandomForestModel`), full demo
  workflows, `compareStrelka*Vcfs.bash` + Q-score deltas quantified

Full results: see `report.md`.

## Tool chains

| | py2 (reference) | py3 (under test) |
|---|---|---|
| Python | CPython 2.7.18 (source build, `--enable-shared --enable-unicode=ucs4`) | CPython 3.13.5 (system) |
| prefix / env | `/workspace/.opencode/state/ref-py27` + `/workspace/.opencode/state/ref-py27-packages` | `/workspace/.opencode/state/training-py3-packages` (pip `--target`) |
| numpy / pandas | 1.15.4 / 0.23.4 | 2.2.6 / 2.3.3 |
| scikit-learn | 0.19.2 | 1.9.1 |
| scipy / joblib | 1.1.0 / 0.13.2 | 1.16.3 / 1.6.0 |
| other | bx-python 0.8.7, threadpool 1.3.2, python-dateutil 2.8.2, pytz 2018.9, six 1.11.0 | threadpoolctl 3.7.0, narwhals 2.27.1, python-dateutil 2.9.0.post0, pytz 2026.5, tzdata 2026.5, six 1.17.0 |
| code | HEAD (unconverted) in `/workspace/repos/strelka` | converted tree in this worktree |
| invocation | `env LD_LIBRARY_PATH=/workspace/.opencode/state/ref-py27/lib PYTHONPATH=/workspace/.opencode/state/ref-py27-packages /workspace/.opencode/state/ref-py27/bin/python2 <script>` | `env PYTHONPATH=/workspace/.opencode/state/training-py3-packages python3 <script>` |

The py2 pins are the wheels of the 2018 era (cp27mu, i.e. UCS-4 unicode —
mandatory: the source build must use `--enable-unicode=ucs4` or the
cp27mu wheels fail with `PyUnicodeUCS4_DecodeUTF8` errors). Install with
`pip install --no-deps` from the wheel list in `design.md` section
"Toolchain py2 de référence (gelée)".

The py3 pins are frozen in
`src/python/scoringModelTraining/requirements-training.txt`
(pandas stays on 2.3.3: `fillna(inplace=…)` is removed in pandas 3.0).

## Data set (`data/`)

Built by `tools/build_parity_dataset.py` (stdlib-only python3,
deterministic — a second build is byte-identical, verified by diff).
Sources, all versioned or derived from versioned files:

* somatic truth: `src/demo/expectedResults/` (17 SNV + 2 indels)
* somatic/germline real calls + features: the demo analyses run with
  `--reportEVSFeatures` (`build/strelkaSomaticDemoAnalysisEVSF/`,
  `build/strelkaGermlineDemoAnalysisEVSF/` — germline single-sample S1)
* synthetic false positives: seeded (`SEED = 20261010`), EVSF values drawn
  data-driven from the observed real feature values (integers: resampled
  from the distinct observed values ±1; floats: observed value ×
  U(0.4, 1.6), zeros kept strict)

Files:

| file | content |
|---|---|
| `somatic_query_snv.vcf` | 17 real TP calls + 50 synthetic FP (pos 100+96i, REF from fasta, ALT = ACGT rotation, QSS 3-40, QUAL 2-12) |
| `somatic_truth_snv.vcf` | the 17 truth calls |
| `somatic_query_indel.vcf` | 2 real TP indels + 20 synthetic FP (pos 150+240i, insertion/deletion alternating) |
| `somatic_truth_indel.vcf` | the 2 truth indels |
| `somatic_fp_regions.bed` | `demo20 8 5008` (labels the synthetic FPs as FP) |
| `germline_query_gvcf.vcf` | mono-sample gVCF of NA12891 (S1), 52 rows: 13 real calls (11 SNP + 2 indels) + 50 synthetic FP (40 SNV + 10 indel, pos 200+90i) |
| `germline_annotated.vcf` | 2-sample annotated VCF (BD:BVT, samples = truth/query) as consumed by `parseAnnotatedTrainingVcf`: 66 rows = 13 TP + 50 FP + 3 FN (TP idx % 5 == 0, POS+1, truth-only, features = empty) |

## Training chain (`tools/run_chain.sh`)

One script, run once per tool chain with the identical data set and
parameter files, so the only difference is the interpreter + the tree:

    bash tools/run_chain.sh <python> <scoringModelTraining-dir> <output-dir>
    # py2: run_chain.sh /workspace/.opencode/state/ref-py27/bin/python2 /workspace/repos/strelka/src/python/scoringModelTraining results/py2
    # py3: run_chain.sh python3 <worktree>/src/python/scoringModelTraining results/py3

Steps (per tool chain):

1. `vcf_to_feature_csv.py` somatic SNV (`--truth`, `--fp-regions`, `--features somatic.snv`)
2. `vcf_to_feature_csv.py` somatic indel
3. `tools/postprocess_indel_csv.py` — adds `NT` (=ref) and `INDELTYPE`
   (1 insertion / 2 deletion) columns consumed by `StrelkaRFIndel.train`
   (not emitted by `vcf_to_feature_csv`); deterministic stdlib script,
   shared by both tool chains (parity is checked on the upstream CSV)
4. `evs_learn.py -m somatic.rf --features somatic.snv -p params/params_somatic_snv.json`
5. `evs_exportmodel.py -c … -v SNV -t 15` → `model_somatic_snv.json`
6. `tools/train_strelka_rf_indel.py` — same flow as `evs_learn` for the
   indel model (shared driver: `strelka_rf_indel` is not imported by
   `evs/__init__` in the somatic tree — pre-existing at HEAD, identical
   under py2 and py3, so `evs_learn.py -m strelka.rf.indel` is not
   selectable by either tool chain)
7. `evs_exportmodel.py -c … -v INDEL` — expected to FAIL on both tool
   chains (pre-existing: `StrelkaRFIndel` has no
   `save_json_strelka_format` + the pickle is not loadable via
   `createFromFile` because the model name is not registered); exit code
   recorded, non-fatal
8. `filterTrainingVcf.py` (germline gVCF → filtered)
9. `parseAnnotatedTrainingVcf.py` (annotated 2-sample VCF → SNV + indel CSV)
10. `evs_learn.py -m germline.rf --features germline.snv -p params/params_germline.json`
11. `evs_learn.py -m germline.rf --features germline.indel -p params/params_germline.json`
12. `evs_exportmodel.py -c … -v SNV -t 15` → `model_germline_snv.json`
13. `evs_exportmodel.py -c … -v INDEL -t 15` → `model_germline_indel.json`

Training parameters (`params/`, identical for both tool chains —
`random_state=42` pinned because it is not set by default; no
`--sample-input` / `--balance_*` so the unseeded `random.sample` of
`evs_learn` is never reached):

* `params_somatic_snv.json` — `{n_jobs: 8, max_depth: 6, n_estimators: 100, random_state: 42}`
* `params_somatic_indel.json` — `{n_jobs: 8, max_depth: 4, n_estimators: 50, max_features: "log2", random_state: 42}`
* `params_germline.json` — `{n_jobs: 8, max_depth: 12, min_samples_leaf: 2, n_estimators: 50, max_features: null, verbose: 10, random_state: 42}`
  (`min_samples_leaf=2` instead of the default 50: the default is
  degenerate on a 54-row data set — documented adaptation, same value on
  both tool chains)

## CP3 — C++ scoring runs

Model JSONs loaded via the workflow options
`--snvScoringModelFile` / `--indelScoringModelFile`
(segment options `--{somatic-}snv-scoring-model-file` /
`--{somatic-}indel-scoring-model-file`), demo data set
(`install/share/demo/strelka/data/`), `runWorkflow.py -m local -j 1 -g 4`:

* `strelkaSomaticCP3m2` / `strelkaSomaticCP3m3` — somatic,
  `--snvScoringModelFile` = JSON_2 / JSON_3 (indel model = shipped in
  both: the somatic indel model is not exportable in either tool chain,
  see step 7)
* `strelkaGermlineCP3m2` / `strelkaGermlineCP3m3` — germline,
  `--snvScoringModelFile` + `--indelScoringModelFile` = JSON_2 / JSON_3
  **without `--exome`**: with `--exome` the germline workflow sets
  `isEVS=False` unconditionally (pre-existing,
  `strelkaGermlineWorkflow.py:412-413`), so the custom models would not be
  applied — genome mode keeps EVS enabled and validates the models

Comparison: `scratch/util/compareStrelka{Somatic,Germline}Vcfs.bash`
(diffs under `results/`) + per-variant Q-score extraction (see report).

## Reproducing

```bash
cd <worktree>
# CP1 + CP2 chain (py2 reference):
env LD_LIBRARY_PATH=/workspace/.opencode/state/ref-py27/lib \
    PYTHONPATH=/workspace/.opencode/state/ref-py27-packages \
    bash scratch/parity/tools/run_chain.sh /workspace/.opencode/state/ref-py27/bin/python2 \
    /workspace/repos/strelka/src/python/scoringModelTraining scratch/parity/results/py2
# CP1 + CP2 chain (py3):
env PYTHONPATH=/workspace/.opencode/state/training-py3-packages \
    bash scratch/parity/tools/run_chain.sh python3 \
    src/python/scoringModelTraining scratch/parity/results/py3
# CP2 analysis:
python3 scratch/parity/tools/compare_models.py \
  py2:py3:scratch/parity/results/py2/model_somatic_snv.json:scratch/parity/results/py3/model_somatic_snv.json \
  py2:py3:scratch/parity/results/py2/model_germline_snv.json:scratch/parity/results/py3/model_germline_snv.json \
  py2:py3:scratch/parity/results/py2/model_germline_indel.json:scratch/parity/results/py3/model_germline_indel.json
# CP2b per-variant scores (4 runs, one per model × tool chain — see report)
# CP3: see above (C++ runs)
# checksums:
sha256sum -c scratch/parity/checksums.sha256
```
