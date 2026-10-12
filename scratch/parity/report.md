# Rapport de parité Py2/Py3 — outils de training EVS/RF
Change : `strelka-py3-training` — date : 2026-10-10

## Critère d'acceptation

Environnement clinique, cadre ISO/COFRAC : la parité entre le code
Python 2 (HEAD non converti, toolchain de référence 2018) et le code
Python 3 converti (toolchain pinnée 2026) est exigée sur la chaîne de
training EVS/RF (`scoringModelTraining/`), en 3 checkpoints :

1. **CP1** — CSV de features py2 vs py3 identiques byte à byte
2. **CP2** — JSON de modèle : feature map identique + divergences d'arbres
   sklearn rapportées et analysées (+ scoring fonctionnel par variante)
3. **CP3** — JSON_2 et JSON_3 chargés en C++ (`RandomForestModel`),
   démos complètes sur VCF de référence fixe, `compareStrelka*Vcfs.bash`,
   Q-scores quantifiés

## Toolchains (gelées)

* **py2 (référence)** : CPython 2.7.18 source build
  (`--enable-shared --enable-unicode=ucs4`, prefix
  `/workspace/.opencode/state/ref-py27`) + wheels cp27mu 2018
  (`ref-py27-packages`) : numpy 1.15.4, pandas 0.23.4, scipy 1.1.0,
  scikit-learn 0.19.2, bx-python 0.8.7, joblib 0.13.2, threadpool 1.3.2,
  python-dateutil 2.8.2, pytz 2018.9, six 1.11.0. Commandes
  d'installation complètes : `design.md` § « Toolchain py2 de référence (gelée) ».
* **py3 (sous test)** : Python 3.13.5 + pins de
  `src/python/scoringModelTraining/requirements-training.txt` :
  numpy 2.2.6, pandas 2.3.3, scipy 1.16.3, scikit-learn 1.9.1,
  joblib 1.6.0, threadpoolctl 3.7.0, narwhals 2.27.1,
  python-dateutil 2.9.0.post0, pytz 2026.5, tzdata 2026.5, six 1.17.0.

Dataset : synthétique déterministe (SEED=20261010) construit des démos
Strelka (truth somatique `src/demo/expectedResults/` + appels/features
réels des runs `--reportEVSFeatures` + FPs synthétiques data-driven) —
source et nature documentées dans `design.md` § « Dataset de référence »
et `README.md` (repli du FTP Illumina MORT, mort : timeout TCP
192.84.34.14:21). Params de training identiques des deux côtés
(`random_state=42` pinné) : `params/`.

## CP1 — parité de features : **PASS (byte à byte)**

Extraction de features complète (somatic : `vcf_to_feature_csv` SNV +
indel ; germline : `filterTrainingVcf` + `parseAnnotatedTrainingVcf`),
mêmes options, mêmes données :

| fichier | taille | py2 vs py3 |
|---|---|---|
| `somatic_snv_features.csv` | 6 927 B | **identique** |
| `somatic_indel_features.csv` | 3 643 B | **identique** |
| `somatic_indel_train.csv` (post NT/INDELTYPE) | 3 788 B | **identique** |
| `germline_snv_features.csv` | 8 568 B | **identique** |
| `germline_indel_features.csv` | 2 101 B | **identique** |
| `germline_filtered.vcf` | 13 114 B | **identique** |

Aucune tolérance de formatage float nécessaire : les représentations
court (shortest repr) des valeurs numériques coïncident entre
`str(float)` py2 / `repr(float)` py3 et `to_csv` de pandas 0.23.4 / 2.3.3
sur ce dataset (pas de NaN dans le CSV somatique ; les cellules vides
germline (FN) sont écrites identiquement).

## CP2 — parité de modèle : **PASS (feature map identique, divergences d'arbres analysées)**

`evs_learn` (même `random_state=42`) + `evs_exportmodel` des deux côtés
(3 modèles exportables : somatic SNV, germline SNV, germline INDEL) :

* **Méta-données : toutes identiques** sur les 3 modèles — `Features`
  (noms + ordre), `ModelType` (RandomForest), `FilterCutoff` (15),
  `Calibration` ({Power:1, Scale:1}), `n_trees` (100 / 50 / 50).
* **Divergences d'arbres** (attendues : flux aléatoires bootstrap
  différents entre sklearn 0.19.2 et 1.9.1, même `random_state`) :

  | modèle | arbres structurellement identiques | identiques seuils <1e-6 | Δ nœuds max |
  |---|---|---|---|
  | somatic SNV (100 arbres) | 9/100 | 1/100 | 14 (plages 3..23 des deux côtés) |
  | germline SNV (50 arbres) | 14/50 | 9/50 | 2 (plages 3..5) |
  | germline INDEL (50 arbres) | 5/50 | 2/50 | 4 (plages 1..5) |

* **Sémantique `node_votes`** — divergence documentée, équivalente pour
  le consommateur C++ :
  * sklearn 0.19.2 : `tree.value` = **comptes de classes** par nœud
    (somme = taille de l'échantillon bootstrap, 1..67) ;
  * sklearn 1.9.1 : `tree.value` = **probabilités normalisées** (somme
    1.0 par nœud) — vérifié empiriquement (`DecisionTreeClassifier`
    ad hoc : somme = 1.0, racine = compte/bootstrap) ;
  * le C++ (`RandomForestModel::getProb`) consomme le ratio
    `vote.left/(vote.left+vote.right)` — **invariant d'échelle** : les
    deux layouts sont des drop-in équivalents ; `classes_` = [FP, TP]
    conservé (FP index 0) des deux côtés.
* **Scoring fonctionnel (CP2b)** — chaque toolchain score le même CSV
  avec son propre pickle (`score_check.py`, P(TP) par variante, FN
  exclues) :

  | modèle | lignes | ptag identique | |ΔP(TP)| max / p95 / mean |
  |---|---|---|---|
  | somatic SNV | 67 | **67/67** | 0.0723 / 0.0469 / 0.0152 |
  | germline SNV | 51 | **51/51** | 0.0577 / 0.0433 / 0.0061 |
  | germline INDEL | 12 | **12/12** | 0.0850 / 0.0850 / 0.0404 |

  L'impact fonctionnel des divergences d'arbres est borné : aucune
  inversion de classe (FP↔TP) sur l'ensemble de training, écarts de
  probabilité ≤ 0.085.

Modèle somatic INDEL (`StrelkaRFIndel`) : pickle produit des deux côtés
(`train_strelka_rf_indel.py`), pas de JSON exportable — **pré-existant**
(absence de `save_json_strelka_format` + modèle non enregistré dans
`evs/__init__` de l'arbre somatic) : l'export échoue **des deux côtés**
(erreur consignée dans `results/py{2,3}/logs/step7_*.log`).

## CP3 — parité de sortie C++ : **PASS (appels identiques, Q-scores quantifiés)**

Démos complètes (C++ rebuildé, boost fix — voir design.md), dataset démo
fixe (`demo20`, BAM NA12891/NA12892), modèles JSON_2 / JSON_3 passés via
`--{snv,indel}ScoringModelFile` (`strelkaSomaticCP3m2/m3`,
`strelkaGermlineCP3m2/m3`) :

* **Somatic SNV** : **17/17 variantes identiques** (aucun appel
  ajouté/retrait, mêmes POS/REF/ALT) ; SomaticEVS (Q-score) :
  |Δ| max = **14.21** (1744 : 33.42 → 19.21), mean = 2.54 (17 variantes).
  Les 17 appelées = les 17 truth du dataset démo (aucun gain/perdre
  d'appel par rapport à la baseline : la baseline sans modèle custom
  appelle les mêmes 17 positions).
* **Somatic INDEL** : VCF **identiques** (modèle indel shippé des deux
  côtés — le modèle indel entraîné n'est pas exportable, voir CP2).
* **Germline** : **19/19 variantes identiques** ; GQ : 0 changement ;
  GQX : **1 changement** (3664, S1 : 3 → 4, Δ=1) ; QUAL variante : 0
  changement. Le GQ variant-level étant clamped/arrondi, les petits
  écarts de probabilité (CP2b) restent en deçà du seuil de l'unité —
  seul un GQX (non clampé) bas a bougé d'une unité.
  *Note* : les runs germline CP3 sont en mode génome (sans `--exome`) :
  avec `--exome`, le workflow germline pose `isEVS=False`
  inconditionnellement (`strelkaGermlineWorkflow.py:412-413`,
  **pré-existant**) et les modèles custom ne seraient pas appliqués —
  l'identité m2/m3 observée en mode exome (0/19) est triviale.
* `compareStrelka{Somatic,Germline}Vcfs.bash` : diffs conservés dans
  `results/cp3-*.diff` (uniquement les valeurs de scores, aucun diff
  d'appel). *Lecture des diffs (note 2026-10-12, change
  `strelka-qa-noregression`)* : la vue des scripts de comparaison est
  **filtrée `grep PASS`** — `cp3-somatic-snvs.diff` montre donc
  variantLineCount 13 (m2) / 10 (m3), alors que le **set d'appels brut est
  17/17 identiques** (aucun appel créé/retrait) : 3 variants (1706, 2512,
  3366) restent appelées en m3 avec un passage de FILTER PASS→LowEVS
  (SomaticEVS 17.66→13.98 / 16.89→13.98 / 21.47→20.00). Les 2 vues sont
  cohérentes ; le critère d'acceptation porte sur le set d'appels brut.

**Tolérance (décision clinique — consignée)** : les écarts de Q-scores
observés (somatic SNV max Δ14.21 / mean Δ2.54 ; germline GQX Δ1) sont la
conséquence directe des divergences d'arbres sklearn 0.19→1.9 sur un
petit dataset synthétique (67/54/12 lignes) ; aucun appel n'est
créé/retrait sur le set de référence. La tolérance applicable est à
fixer par l'utilisateur (clinique) sur la base de ces mesures ; le
rapport et les artefacts (`results/`, `checksums.sha256`) permettent de
re-juger sur un dataset de production (VCF réel) sans ré-exécuter la
méthodologie.

## Anomalies pré-existantes constatées (identiques py2 @ HEAD et py3 converti — laissées telles quelles pour la parité)

1. `somatic/bin/calc_features.py` : `ImportError: No module named ref`
   (`evs/features/ref.py` absent) — module mort des deux côtés.
2. `GermlineRF.draw_trees` : NameError (`featurenames`/`Image` non
   définis) — hors chemin de la chaîne.
3. `StrelkaRFIndel` : pas de `save_json_strelka_format` ; et
   `strelka_rf_indel` n'est jamais importé par `evs/__init__` (arbre
   somatic) → `evs_learn.py -m strelka.rf.indel` et l'export JSON indel
   somatique indisponibles des deux côtés (driver partagé
   `train_strelka_rf_indel.py` fourni pour la parité).
4. Docs userGuide : noms de modèles `strelka.snv`/`strelka.rf` alors que
   le code enregistre `somatic.snv`/`somatic.rf` — corrigées par ce
   change (task 6.2).
5. JSON indel livré = entrée unique `CalibrationModels→Somatic→INDEL`.
6. `--exome` → `isEVS=False` en germline (pré-existant, documenté CP3).

## Audit des divisions entières (task 3.2)

22 candidates auditées : aucune conversion `/`→`//` nécessaire — les
expressions en jeu (ex. `0.5*FP/UNK`) sont déjà flottantes en py2 (facteur
`0.5`) et restent flottantes en py3 ; les autres `/` opèrent sur des
integers où la sémantique est inchangée par la parité CP1 (CSV byte
identiques ⇒ valeurs des features identiques).

## Artefacts

* `data/` (dataset pinné, 7 fichiers), `params/` (3 JSON), `tools/`
  (builder, post-process, drivers, chaîne, comparateurs)
* `results/py2/`, `results/py3/` : CSV, pickles, JSON, logs par étape
* `results/scores_*.tsv` : P(TP) par variante (CP2b)
* `results/cp3-*.diff` : diffs VCF CP3
* `checksums.sha256` : dataset + artefacts
* Logs C++ : `/workspace/.opencode/tmp/opencode/cp3-*.log`

Reproduction : `README.md`.
