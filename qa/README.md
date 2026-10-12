# Strelka — dispositif qualité de non-régression (`qa/`)

Changement : `strelka-qa-noregression` (store OpenSpec `/workspace/openspec`).

Objectif : dans un environnement clinique (cadre **EN ISO/IEC 15189:2022** +
**EN ISO/IEC 27001:2022**, accréditation COFRAC — adjacentes selon périmètre :
ISO/IEC 17025:2017, IEC 62304/MDR 2017/745), chaque modification de l'outil
(code, dépendances, toolchain, données) doit démontrer l'absence de régression
de qualité — mêmes variants appelés, Q-scores dans les tolérances — avec des
preuves **cumulatives, horodatées, hachées et traçables** (« au fur et à
mesure »).

**Séparation des rôles** : le repo tient le *dossier technique* (preuves
reproductibles, hachées, re-jeuables — le re-run de la porte re-joue la
démo). Le *QMS du labo* tient les procédures/approvals et y pointe
(`qa/evidence/<rapport>` + changes archivés OpenSpec). Le repo ne contient
pas de contenu QMS.

## Baseline golden (consolidée, référencée en place — pas dupliquée)

| artefact | chemin (repo) | ce qu'il protège | checksum |
|---|---|---|---|
| dataset de parité (7 fichiers) | `scratch/parity/data/` | parité training Py2/Py3 (CP1/CP2) | `scratch/parity/checksums.sha256` (re-jeu par la porte) |
| résultats attendus démos somatique (caller) | `src/demo/expectedResults/` | qualité des appels somatique (la démo somatique compare elle-même en step 3) | hash agrégé dans chaque rapport |
| résultats attendus démos germline (caller) | `src/demo/expectedResultsGermline/` (4 VCF, métadonnées variables filtrées) | qualité des appels germline — la comparaison est faite **par la porte** (`demo_germline_golden`) : la step 3 de comparaison du script de démo germline est commentée en amont (pré-existant, aucune baseline germline n'existait) | hash agrégé dans chaque rapport |
| résultats attendus démos mito (caller) | `src/demo/expectedResultsMito/` | qualité des appels mito | versionné (git) |
| lockfile training py3 | `src/python/scoringModelTraining/requirements-training.txt` | environnement de training py3 pinné | sha256 dans chaque rapport |
| toolchain C++ pinnée | `redist/` (cmake 3.31.12, Boost 1.92.0 subset, htslib/samtools 1.24, …) | build C++ reproductible | hash agrégé dans chaque rapport |
| rapport de parité (baseline clinique) | `scratch/parity/report.md` | état de référence Py2/Py3 validé 2026-10-10 | versionné (git) |

Déterminisme de la baseline germline vérifié à la création (2026-10-12) :
2 runs du démo germline → VCF identiques byte-à-byte après filtrage des
métadonnées variables (`fileDate`, `source_version`, `startTime`,
`reference`, `cmdline` — filtre identique au `filterVariableMetadata` du
script de démo). Toute divergence future = la porte FAIL (`demo_germline_golden`).

Toute modification d'un artefact de baseline ⇒ son checksum change ⇒ la porte
le signale (étape `baseline_checksums` FAIL) ⇒ le diff entre deux rapports de
`qa/evidence/` montre exactement l'input modifié (c'est la dérive, visible
d'un coup d'œil).

## Inventaire des toolchains (27001 A.5.9)

### py3 — environnement de training (pinné)
`src/python/scoringModelTraining/requirements-training.txt` (numpy 2.2.6,
pandas 2.3.3, scipy 1.16.3, scikit-learn 1.9.1, joblib 1.6.0,
threadpoolctl 3.7.0, narwhals 2.27.1, python-dateutil 2.9.0.post0,
pytz 2026.5, tzdata 2026.5, six 1.17.0 — pandas reste 2.3.3 :
`fillna(inplace=…)` retirée en pandas 3.0). Installation :
`pip install --target <dir> -r requirements-training.txt` (container :
`/workspace/.opencode/state/training-py3-packages`, variable
`STRELKA_QA_PY3_PKGS`).

### py2 — toolchain de référence (gelée, re-constituible from-scratch)
État du conteneur de build, **pas** du code versionné ; re-jeu vérifié le
2026-10-10 (source : change archivé `2026-10-10-strelka-py3-training`,
section « Toolchain py2 de référence (gelée) ») :

1. CPython **2.7.18** build source, tarball officiel python.org
   (`/workspace/.opencode/state/ref/Python-2.7.18.tgz`) :
   `./configure --prefix=<prefix> --enable-shared --enable-unicode=ucs4`
   (`--enable-unicode=ucs4` **obligatoire** — les wheels cp27mu sont UCS-4)
   + `make -j16 && make install`
2. `python2 -m ensurepip`
3. pins époque 2018 (contemporains de Strelka 2.9.10, nov. 2018),
   `pip install --target <pkgs> --no-deps` :

| wheel (PyPI, 2026-10-10) | sha256 |
|---|---|
| numpy-1.15.4-cp27-cp27mu-manylinux1_x86_64.whl | 0df89ca13c25eaa1621a3f09af4c8ba20da849692dcae184cb55e80952c453fb |
| pandas-0.23.4-cp27-cp27mu-manylinux1_x86_64.whl | be4715c9d8367e51dbe6bc6d05e205b1ae234f0dc5465931014aa1c4af44c1ba |
| scipy-1.1.0-cp27-cp27mu-manylinux1_x86_64.whl | 08237eda23fd8e4e54838258b124f1cd141379a5f281b0a234ca99b38918c07a |
| scikit_learn-0.19.2-cp27-cp27mu-manylinux1_x86_64.whl | a402c1484fe65df42d5dbc22a58e0695fe3afe2b0b229aee2a09c6d60ba8e5c2 |
| bx_python-0.8.7-cp27-cp27mu-manylinux2010_x86_64.whl | 21a58a9da6da1ae1ce1171aae55f990f1b3afaa17ecee4dfcb3d3dd51252c617 |
| joblib-0.13.2-py2.py3-none-any.whl | 21e0c34a69ad7fde4f2b1f3402290e9ec46f545f15f1541c582edfe05d87b63a |
| threadpool-1.3.2-py2-none-any.whl | 826e0c25d4cd3c1e58c2e16a6219f60ba13b5b34a5837a7f8659712f88cb7cdd |
| python_dateutil-2.8.2-py2.py3-none-any.whl | 961d03dc3453ebbc59dbdea9e4e11c5651520a876d0f4db161e8674aae935da9 |
| pytz-2018.9-py2.py3-none-any.whl | 32b0891edff07e28efe91284ed9c31e123d84bea3fd98e1f72be2508f43ef8d9 |
| six-1.11.0-py2.py3-none-any.whl | 832dc0e10feb1aa2c68dcc57dbb658f1c7e65b9b61af69048abc87a2db00a0eb |

   (pandas 0.23.4 : wheel local de repli conservé sous
   `/workspace/.opencode/state/ref/pandas-0.23.4-cp27mu.whl` si PyPI devient
   indisponible.)
4. Invocation : `env LD_LIBRARY_PATH=<prefix>/lib PYTHONPATH=<pkgs> <prefix>/bin/python2 <script>`
   (build `--enable-shared` ⇒ `LD_LIBRARY_PATH` requis au runtime).

Localisations container (défauts de `qa/validate.sh`, surchargeables par
variables d'env) : prefix `/workspace/.opencode/state/ref-py27`, packages
`/workspace/.opencode/state/ref-py27-packages`.

### py2 — arbre de référence non converti
Le code non converti (état `python3` pré-merge du change
`strelka-py3-training`) est **matérialisé à la volée** par la porte :
`git archive 2cf5f0fc src/python/scoringModelTraining` (SHA pinné =
`STRELKA_QA_PY2_TREE_SHA`, consigné dans chaque rapport). Aucune copie hors
git, aucune dérive possible.

### C++
gcc (container Debian 13 : 14.2.0), cmake **3.31.12** + **Boost 1.92.0**
(subset) + htslib/samtools **1.24** : tarballs pinnés dans `redist/`,
bootstrap par `./configure` (build séparé : `mkdir build && cd build &&
../configure --prefix=<install> --build-type=Release --jobs=16 && make -j16 &&
make install`). Versions effectives consignées dans chaque rapport
(`gcc`, `cmake`, `boost`).

## La porte (`qa/validate.sh`)

Une commande ; chaque run produit un rapport dans `qa/evidence/` :

```
qa/validate.sh            # fast
qa/validate.sh --full     # + parité Py2/Py3 CP1/CP2/CP2b
qa/validate.sh --cp3      # + runs C++ de scoring avec les modèles parité
qa/validate.sh --install DIR --keep-work
```

| couche | étapes | critère de FAIL |
|---|---|---|
| fast | `compileall_training` (module `scoringModelTraining`) · `unittest_training` (suite `tests/`) · `baseline_checksums` (re-jeu `scratch/parity/checksums.sha256`) · `demo_somatic` (script installé vs `expectedResults`, comparaison native step 3 du script) · `demo_germline` (script installé) + `demo_germline_golden` (comparaison par la porte vs `expectedResultsGermline`, métadonnées variables filtrées) | rc ≠ 0 de n'importe quelle étape |
| full | `parity_ref_tree` (matérialisation `git archive` du SHA pinné) · `parity_py2_chain` + `parity_py3_chain` (`scratch/parity/tools/run_chain.sh`, dataset pinné) · `parity_cp1` (6 CSV/VCF **byte-identiques** — critère) · `parity_cp2` (analyse arbres sklearn, `compare_models.py` — quantifié, informationnel) · `parity_cp2b` (P(TP) par variante avec le modèle de chaque toolchain, `score_check.py` — quantifié, informationnel) | FAIL si une chaîne échoue ou CP1 diverge ; CP2/CP2b rapportent les écarts (les divergences d'arbres sklearn 0.19→1.9 sont attendues et documentées dans `scratch/parity/report.md` — la décision clinique porte sur ces mesures) |
| cp3 | `cp3_somatic_m2/m3` + `cp3_germline_m2/m3` (runs C++ complets, modèles JSON_2/JSON_3 via `--{snv,indel}ScoringModelFile`, démo20 ; germline en **mode génome** — `--exome` force `isEVS=False` pré-existant) · `cp3_check` (**variant-set identique** = critère ; deltas de scores quantifiés) | FAIL si une variante est créée/retraitée ; les deltas Q-scores sont consignés dans la ligne de résultat |

**N/A explicite** : une couche dont la précondition manque (toolchain py2
absente, install tree absent…) n'est jamais silencieuse — le rapport porte
`N/A (explicit)` + la raison (ex. « py2 reference tool chain absent »).

Variables d'env : `STRELKA_QA_PY2_BIN`, `STRELKA_QA_PY2_LIB`,
`STRELKA_QA_PY2_PKGS`, `STRELKA_QA_PY3_PKGS`, `STRELKA_QA_PY2_TREE_SHA`
(défauts = localisations container, voir haut).

Coût d'exécution (container, 16 cœurs) : fast ≈ 10 min (2 démos) ; `--full`
≈ +15 min (chaînes ×2) ; `--cp3` ≈ +15 min (4 runs C++).

## Registre de preuves (`qa/evidence/`)

- **Append-only** : un run = `<UTC>-<gitsha7>.{md,json}` (suffixe `-n` en cas
  de collision de timestamp) ; aucun rapport existant n'est modifié.
- **Auto-haché** : le `.md` porte le sha256 du `.json` — vérifiable :
  `sha256sum -c` sur le `.json` comparé à la valeur du `.md`.
- **Auto-portant** : chaque rapport consigne le niveau de porte exécuté, le
  résultat et les *key numbers* de chaque étape (ligne `result`), l'inventaire
  d'environnement (git, python3, gcc, cmake, boost, toolchains py2/py3, SHA de
  l'arbre de référence) et les sha256 des inputs pinnés (baseline, lockfile
  py3, expectedResults, redist). Les logs complets restent dans le work
  temporaire (`--keep-work` pour les conserver) ; leur sha256 est dans le
  rapport.
- **Dérive** : `diff` entre deux rapports = uniquement les inputs réellement
  changés entre les deux runs (horodatage, résultat, hachés).

## Fiche qualité par change (template)

À compléter dans le `design.md` de chaque change OpenSpec qui touche l'outil
:

```markdown
### Qualité & non-régression
- Périmètre (fichiers/surfaces touchés) :
- Risques identifiés (scoring, appels, pins, données) :
- Preuves : qa/evidence/<UTC>-<gitsha7>.md (niveau fast/full/cp3)
  + rapports spécifiques (ex. scratch/parity/report.md)
- Résidus acceptés & tolérances cliniques (décision, auteur, date) :
- Approbation (qui, quand) :
```

Un change archivé sans preuve de non-régression doit le justifier dans sa
fiche (ex. changement purement documentaire). Rattrapage des changes
mergés le 2026-10-09/10 : `strelka-py3-training` conforme (rapport parity
commité) ; `strelka-py3-runtime` / `strelka-redist-modern` validés par démos
E2E + ctest (à pointer dans leur fiche au prochain rattrapage).

## Triggers et périodicité

- **Trigger** : tout change code / dépendances / toolchain / dataset ⇒
  porte **fast** au minimum **avant merge** ; `--full` quand la surface
  touche les outils de training ou les pins ; `--cp3` quand la surface
  touche le scoring C++ ou les modèles (validation de parité).
- **Périodique** (15189 §8.2 revue de la direction / 17025 7.8 revue
  périodique) : re-run annuel de la porte (fast au minimum).
- **Health check toolchain** : réinstallation from-scratch de la référence
  py2 (procédure ci-dessus, effectuée 2026-10-10 — re-vérifier à chaque
  changement majeur du container) + vérification que les pins py3
  s'installent dans un environnement neuf.
- Tout **changement de baseline** (dataset, expectedResults, pins, redist)
  ⇒ `--full`/`--cp3` + mise à jour des checksums + entrée dans le registre.

## Mapping normes

| norme / clause | ce que le dispositif fournit |
|---|---|
| ISO 15189 §7.3 (validation du processus) | porte `qa/validate.sh` + niveaux fast/full/cp3 consignés |
| ISO 15189 §7.5 (enregistrements traçables) | `qa/evidence/` append-only auto-haché, inputs + environnement hachés |
| ISO 15189 §8.2 (revue périodique) | périodicité annuelle + health check toolchain (ci-dessus) |
| ISO 27001 A.5.9 (inventaire des actifs) | inventaire des toolchains pinnées + checksums (ci-dessus) |
| ISO 27001 A.8.25 (intégrité) | sha256 des inputs pinnés + du rapport lui-même ; dérive = diff visible |
| ISO 27001 A.8.30 (gestion des changements) | trigger « tout change ⇒ porte avant merge » + fiche qualité par change |
| ISO 17025 7.8 (revue périodique) | idem 15189 §8.2 (accréditation COFRAC) |
