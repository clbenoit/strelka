#!/bin/bash
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT OTHER WARRANTY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

# Runs the full EVS training chain on a fixed data set, for one tool chain.
#
# The same script and the same data set are used for the py2 and py3
# reference runs, so the only difference between the two is the interpreter
# (and its third-party pins) and the tree under test:
#   py2 : CPython 2.7.18 ucs4 reference tool chain, unconverted code (HEAD)
#   py3 : Python 3.13 + pinned training env, converted code (this worktree)
#
# Usage:
#   run_chain.sh <python-interpreter> <scoringModelTraining-dir> <output-dir>
#
# The somatic indel model export (step 7) is expected to FAIL on both
# tool chains (pre-existing: StrelkaRFIndel has no
# save_json_strelka_format); the exit code is recorded, not fatal.

set -u
set -o pipefail

if [ $# -ne 3 ]; then
    echo "usage: $0 <python-interpreter> <scoringModelTraining-dir> <output-dir>" >&2
    exit 2
fi

PY="$1"
TREE="$2"
OUT="$3"

HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="$HERE/../data"
PARAMS="$HERE/../params"
TOOLS="$HERE"

mkdir -p "$OUT/logs"
SOMBIN="$TREE/somatic/bin"
GERMBIN="$TREE/germline/bin"

step=0
run() {
    # run <step-name> <cmd...>
    local name="$1"; shift
    step=$((step + 1))
    echo "==== [$step] $name"
    echo "+ $*"
    if "$@" 2>&1 | tee "$OUT/logs/step${step}_$name.log"; then
        return 0
    else
        local rc=$?
        echo "!! step $name exited $rc" >&2
        return $rc
    fi
}

run 1_somatic_snv_csv \
    "$PY" "$SOMBIN/vcf_to_feature_csv.py" "$DATA/somatic_query_snv.vcf" \
        -o "$OUT/somatic_snv_features.csv" \
        --truth "$DATA/somatic_truth_snv.vcf" \
        --fp-regions "$DATA/somatic_fp_regions.bed" \
        --features somatic.snv || exit 1

run 2_somatic_indel_csv \
    "$PY" "$SOMBIN/vcf_to_feature_csv.py" "$DATA/somatic_query_indel.vcf" \
        -o "$OUT/somatic_indel_features.csv" \
        --truth "$DATA/somatic_truth_indel.vcf" \
        --fp-regions "$DATA/somatic_fp_regions.bed" \
        --features somatic.indel || exit 1

run 3_somatic_indel_postprocess \
    "$PY" "$TOOLS/postprocess_indel_csv.py" \
        "$OUT/somatic_indel_features.csv" "$OUT/somatic_indel_train.csv" || exit 1

run 4_somatic_snv_learn \
    "$PY" "$SOMBIN/evs_learn.py" -m somatic.rf --features somatic.snv \
        -p "$PARAMS/params_somatic_snv.json" \
        -o "$OUT/model_somatic_snv.pkl" \
        "$OUT/somatic_snv_features.csv" || exit 1

run 5_somatic_snv_export \
    "$PY" "$SOMBIN/evs_exportmodel.py" -c "$OUT/model_somatic_snv.pkl" \
        -o "$OUT/model_somatic_snv.json" -v SNV -t 15 || exit 1

# strelka_rf_indel is not registered by the somatic evs/__init__ (pre-existing,
# identical under py2 and py3), so evs_learn.py cannot select it; a shared
# driver runs the same flow instead (see tools/train_strelka_rf_indel.py).
run 6_somatic_indel_learn \
    "$PY" "$TOOLS/train_strelka_rf_indel.py" \
        --libdir "$SOMBIN/../lib" \
        -p "$PARAMS/params_somatic_indel.json" \
        -o "$OUT/model_somatic_indel.pkl" \
        "$OUT/somatic_indel_train.csv" || exit 1

# Expected to fail on both tool chains: createFromFile needs "strelka.rf.indel"
# registered (same pre-existing gap) and StrelkaRFIndel has no
# save_json_strelka_format. The exit code is recorded, not fatal.
run 7_somatic_indel_export \
    "$PY" "$SOMBIN/evs_exportmodel.py" -c "$OUT/model_somatic_indel.pkl" \
        -o "$OUT/model_somatic_indel.json" -v INDEL -t 15
echo "somatic indel export rc=$? (expected: non-zero, pre-existing)"

run 8_germline_filter \
    bash -c "$PY '$GERMBIN/filterTrainingVcf.py' < '$DATA/germline_query_gvcf.vcf' > '$OUT/germline_filtered.vcf'" || exit 1

run 9_germline_parse \
    bash -c "$PY '$GERMBIN/parseAnnotatedTrainingVcf.py' < '$DATA/germline_annotated.vcf' --snvOutput '$OUT/germline_snv_features.csv' --indelOutput '$OUT/germline_indel_features.csv'" || exit 1

run 10_germline_snv_learn \
    "$PY" "$GERMBIN/evs_learn.py" -m germline.rf --features germline.snv \
        -p "$PARAMS/params_germline.json" \
        -o "$OUT/model_germline_snv.pkl" \
        "$OUT/germline_snv_features.csv" || exit 1

run 11_germline_indel_learn \
    "$PY" "$GERMBIN/evs_learn.py" -m germline.rf --features germline.indel \
        -p "$PARAMS/params_germline.json" \
        -o "$OUT/model_germline_indel.pkl" \
        "$OUT/germline_indel_features.csv" || exit 1

run 12_germline_snv_export \
    "$PY" "$GERMBIN/evs_exportmodel.py" -c "$OUT/model_germline_snv.pkl" \
        -o "$OUT/model_germline_snv.json" -v SNV -t 15 || exit 1

run 13_germline_indel_export \
    "$PY" "$GERMBIN/evs_exportmodel.py" -c "$OUT/model_germline_indel.pkl" \
        -o "$OUT/model_germline_indel.json" -v INDEL -t 15 || exit 1

echo "chain complete: $OUT"
