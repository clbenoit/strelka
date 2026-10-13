#!/usr/bin/env bash
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
# Strelka non-regression gate (change strelka-qa-noregression).
#
# One command, layered; one self-hashed report per run in qa/evidence/:
#
#   fast (default)  : compileall + unittest + baseline checksums
#                     + caller golden tests (somatic/germline demos vs
#                     expectedResults, via the installed demo scripts)
#   --full          : + Py2/Py3 parity CP1/CP2/CP2b (reference tool chains;
#                     each layer consigné N/A when its tool chain is absent,
#                     never silent)
#   --cp3           : + C++ scoring runs with the parity model JSONs +
#                     semantic check (variant-set parity = pass criterion,
#                     score deltas quantified)
#
#   --install DIR   : Strelka install tree (default: <repo>/install)
#   --keep-work     : keep the temp work dir (full step logs)
#
# Environment overrides (defaults = container reference locations):
#   STRELKA_QA_PY2_BIN    /workspace/.opencode/state/ref-py27/bin/python2
#   STRELKA_QA_PY2_LIB    /workspace/.opencode/state/ref-py27/lib
#   STRELKA_QA_PY2_PKGS   /workspace/.opencode/state/ref-py27-packages
#   STRELKA_QA_PY3_PKGS   /workspace/.opencode/state/training-py3-packages
#   STRELKA_QA_PY2_TREE_SHA  2cf5f0fcdad70550bee660243ab31173509ef98c
#        (python3 @ pre-merge of change 3 = unconverted training tree,
#         materialized with `git archive`; see qa/README.md)
#
# Procedure, triggers and norm mapping: qa/README.md.

set -u
set -o pipefail

QA_DIR=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$QA_DIR/.." && pwd)
QA_TOOLS="$QA_DIR/tools"

FULL=0
CP3=0
KEEP_WORK=0
INSTALL_DIR="$ROOT/install"

usage() {
    cat <<'EOF'
qa/validate.sh [--full] [--cp3] [--install DIR] [--keep-work]

  (no option)  fast: compileall + unittest + baseline checksums
               + caller golden tests (demos vs expectedResults)
  --full       + Py2/Py3 parity CP1/CP2/CP2b (N/A consigné if tool chains
               absent)
  --cp3        + C++ scoring runs with the parity model JSONs + semantic
               check
  --install    Strelka install tree (default: <repo>/install)
  --keep-work  keep the temp work dir (full step logs)
EOF
}
while [ $# -gt 0 ]; do
    case "$1" in
        --full) FULL=1 ;;
        --cp3) CP3=1; FULL=1 ;;
        --keep-work) KEEP_WORK=1 ;;
        --install) INSTALL_DIR=$2; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown argument: $1 (try --help)" >&2; exit 2 ;;
    esac
    shift
done
if [ $CP3 -eq 1 ]; then
    LEVEL=cp3
elif [ $FULL -eq 1 ]; then
    LEVEL=full
else
    LEVEL=fast
fi

PY2_BIN=${STRELKA_QA_PY2_BIN:-/workspace/.opencode/state/ref-py27/bin/python2}
PY2_LIB=${STRELKA_QA_PY2_LIB:-/workspace/.opencode/state/ref-py27/lib}
PY2_PKGS=${STRELKA_QA_PY2_PKGS:-/workspace/.opencode/state/ref-py27-packages}
PY3_PKGS=${STRELKA_QA_PY3_PKGS:-/workspace/.opencode/state/training-py3-packages}
PY2_TREE_SHA=${STRELKA_QA_PY2_TREE_SHA:-2cf5f0fcdad70550bee660243ab31173509ef98c}

TS=$(date -u +%Y%m%dT%H%M%SZ)
GITSHA=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo nogit)
STAMP="${TS}-${GITSHA}"
WORK=$(mktemp -d -t qa-validate.XXXXXX)
STEPS="$WORK/steps.tsv"
: > "$STEPS"

cleanup() {
    if [ $KEEP_WORK -eq 1 ]; then
        echo "work dir kept: $WORK"
    else
        rm -rf "$WORK"
    fi
}
trap cleanup EXIT

# step <name> <result-grep-pattern|-> <cmd...>
# Runs the command, logs to $WORK/<name>.log, appends
# name<TAB>rc<TAB>log<TAB>result to the steps TSV. `result` = first line of
# the log matching the pattern ("" when pattern is "-"), so the report
# carries the key numbers of each step (self-contained evidence).
step() {
    local name=$1 pattern=$2; shift 2
    echo ">> qa: $name"
    local rc=0
    "$@" > "$WORK/$name.log" 2>&1 || rc=$?
    local result=""
    if [ -n "$pattern" ] && [ "$pattern" != "-" ]; then
        result=$(grep -E "$pattern" "$WORK/$name.log" 2>/dev/null \
            | head -1 | tr -d '\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    fi
    printf '%s\t%s\t%s\t%s\n' "$name" "$rc" "$WORK/$name.log" "$result" \
        >> "$STEPS"
    if [ $rc -eq 0 ]; then
        echo "   OK: $name${result:+ — $result}"
    else
        echo "   FAIL($rc): $name (log: $WORK/$name.log)"
    fi
}

# step_na <name> <note>  — explicit skip, never silent
step_na() {
    local name=$1 note=$2
    printf '%s\n' "$note" > "$WORK/$name.log"
    printf '%s\t%s\t%s\t%s\n' "$name" -1 "$WORK/$name.log" "$note" >> "$STEPS"
    echo "   N/A: $name ($note)"
}

has_py2() { [ -x "$PY2_BIN" ] && [ -d "$PY2_PKGS" ]; }
has_py3env() { [ -d "$PY3_PKGS" ]; }
has_install() { [ -x "$INSTALL_DIR/bin/runStrelkaSomaticWorkflowDemo.bash" ] \
    && [ -x "$INSTALL_DIR/bin/runStrelkaGermlineWorkflowDemo.bash" ]; }

PY3_TREE="$ROOT/src/python/scoringModelTraining"
DEMO_DATA=""
if has_install; then
    DEMO_DATA="$INSTALL_DIR/share/demo/strelka/data"
fi

# ---------------------------------------------------------------- fast layer

step compileall_training - \
    python3 -m compileall -q "$PY3_TREE"
# no output on success: the log being empty + rc 0 is the result

step unittest_training '^(Ran|OK|FAILED)' bash -c \
    'cd "$1" && env PYTHONPATH="${2:-}" python3 -m unittest discover -s . -t . -v' \
    _ "$PY3_TREE/tests" "$PY3_PKGS"

step baseline_checksums 'files re-verified' bash -c \
    'cd "$1" && n=$(grep -c "" checksums.sha256) && sha256sum -c checksums.sha256 && echo "$n/$n files re-verified"' \
    _ "$ROOT/scratch/parity"

if has_install; then
    mkdir -p "$WORK/demo-somatic" "$WORK/demo-germline"
    step demo_somatic 'No differences between expected' \
        bash -c \
        'cd "$1" && env PATH="$2/bin:$PATH" "$2/bin/runStrelkaSomaticWorkflowDemo.bash"' \
        _ "$WORK/demo-somatic" "$INSTALL_DIR"
    step demo_germline 'Demo/verification successfully completed' \
        bash -c \
        'cd "$1" && env PATH="$2/bin:$PATH" "$2/bin/runStrelkaGermlineWorkflowDemo.bash"' \
        _ "$WORK/demo-germline" "$INSTALL_DIR"
    # the germline demo script's own comparison step is commented out
    # upstream (pre-existing; no germline expected results existed) — the
    # golden check is performed here against the versioned baseline
    # src/demo/expectedResultsGermline/ (variable metadata filtered, as in
    # the upstream filterVariableMetadata)
    step demo_germline_golden 'no differences' bash -c '
        set -eu
        E="$1"; R="$2"
        filt() { awk "!/^##(fileDate|source_version|startTime|reference|cmdline)/"; }
        n=0
        for f in $(ls "$E"); do
            diff <(gzip -dc "$E/$f" | filt) <(gzip -dc "$R/$f" | filt)
            n=$((n+1))
        done
        echo "$n/$n germline files: no differences (variable metadata filtered)"' \
        _ "$ROOT/src/demo/expectedResultsGermline" \
          "$WORK/demo-germline/strelkaGermlineDemoAnalysis/results/variants"

    # ---------------- python3 compatibility matrix (change
    # strelka-pyflow-py3-threading-fix): the vendored pyflow 1.1.20
    # collides with CPython 3.8-3.12 private names pre-fix (StoppableThread
    # self._stop shadows Thread._stop()); the matrix proves the fix on
    # pinned interpreters, N/A consigné when one is absent.
    step py3_smoke_native 'pyflow threading smoke: OK' bash -c \
        'python3 "$1" "$2"' \
        _ "$ROOT/qa/tools/pyflow_threading_smoke.py" "$INSTALL_DIR"
    for spec in "3.10:PY310_BIN" "3.11:PY311_BIN" "3.12:PY312_BIN"; do
        ver=${spec%%:*}; var=${spec##*:}
        pybin=${!var}
        sname=py3_smoke_${ver//./}
        if [ -n "$pybin" ] && [ -x "$pybin" ]; then
            step $sname 'pyflow threading smoke: OK' bash -c \
                '"$1" "$2" "$3"' _ "$pybin" \
                "$ROOT/qa/tools/pyflow_threading_smoke.py" "$INSTALL_DIR"
        else
            step_na $sname "python $ver interpreter absent (STRELKA_QA_${var})"
        fi
    done
    if [ -n "${PY312_BIN:-}" ] && [ -x "$PY312_BIN" ]; then
        mkdir -p "$WORK/demo312-somatic" "$WORK/demo312-germline"
        step py312_demo_somatic 'No differences between expected' bash -c '
            cd "$1" && env PATH="$(dirname "$3"):$2/bin:$PATH" \
                "$2/bin/runStrelkaSomaticWorkflowDemo.bash"' \
            _ "$WORK/demo312-somatic" "$INSTALL_DIR" "$PY312_BIN"
        step py312_demo_germline 'Demo/verification successfully completed' bash -c '
            cd "$1" && env PATH="$(dirname "$3"):$2/bin:$PATH" \
                "$2/bin/runStrelkaGermlineWorkflowDemo.bash"' \
            _ "$WORK/demo312-germline" "$INSTALL_DIR" "$PY312_BIN"
    else
        step_na py312_demo_somatic "python 3.12 interpreter absent (STRELKA_QA_PY312_BIN)"
        step_na py312_demo_germline "python 3.12 interpreter absent (STRELKA_QA_PY312_BIN)"
    fi
else
    step_na demo_somatic "install tree absent ($INSTALL_DIR)"
    step_na demo_germline "install tree absent ($INSTALL_DIR)"
    step_na demo_germline_golden "install tree absent ($INSTALL_DIR)"
    step_na py3_smoke_native "install tree absent ($INSTALL_DIR)"
    step_na py3_smoke_310 "install tree absent ($INSTALL_DIR)"
    step_na py3_smoke_311 "install tree absent ($INSTALL_DIR)"
    step_na py3_smoke_312 "install tree absent ($INSTALL_DIR)"
    step_na py312_demo_somatic "install tree absent ($INSTALL_DIR)"
    step_na py312_demo_germline "install tree absent ($INSTALL_DIR)"
fi

# ----------------------------------------------------------- full (parity) layer
if [ $FULL -gt 0 ]; then
    step parity_ref_tree '^[0-9]+$' bash -c \
        'set -e; mkdir -p "$3"; git -C "$1" archive "$2" src/python/scoringModelTraining | tar -x -C "$3"; find "$3/src/python/scoringModelTraining" -name "*.py" | wc -l' \
        _ "$ROOT" "$PY2_TREE_SHA" "$WORK/ref"
    REF_TREE="$WORK/ref/src/python/scoringModelTraining"

    if has_py2; then
        step parity_py2_chain 'somatic indel export rc=[0-9]+' \
            env LD_LIBRARY_PATH="$PY2_LIB" PYTHONPATH="$PY2_PKGS" \
            bash "$ROOT/scratch/parity/tools/run_chain.sh" \
            "$PY2_BIN" "$REF_TREE" "$WORK/parity/py2"
    else
        step_na parity_py2_chain "py2 reference tool chain absent ($PY2_BIN / $PY2_PKGS)"
    fi

    if has_py3env; then
        step parity_py3_chain 'somatic indel export rc=[0-9]+' \
            env PYTHONPATH="$PY3_PKGS" \
            bash "$ROOT/scratch/parity/tools/run_chain.sh" \
            python3 "$PY3_TREE" "$WORK/parity/py3"
    else
        step_na parity_py3_chain "py3 pinned env absent ($PY3_PKGS)"
    fi

    if [ -d "$WORK/parity/py2" ] && [ -d "$WORK/parity/py3" ] \
        && [ -f "$WORK/parity/py2/model_germline_indel.json" ] \
        && [ -f "$WORK/parity/py3/model_germline_indel.json" ]; then
        step parity_cp1 'files byte-identical' bash -c '
            set -e
            cd "$1"
            ok=0
            for f in somatic_snv_features.csv somatic_indel_features.csv \
                     somatic_indel_train.csv germline_filtered.vcf \
                     germline_snv_features.csv germline_indel_features.csv; do
                cmp -s "py2/$f" "py3/$f" && ok=$((ok+1)) || echo "DIFF: $f"
            done
            echo "$ok/6 files byte-identical"
            test $ok -eq 6' _ "$WORK/parity"

        step parity_cp2 'structurally identical' \
            python3 "$ROOT/scratch/parity/tools/compare_models.py" \
            "py2:py3:$WORK/parity/py2/model_somatic_snv.json:$WORK/parity/py3/model_somatic_snv.json" \
            "py2:py3:$WORK/parity/py2/model_germline_snv.json:$WORK/parity/py3/model_germline_snv.json" \
            "py2:py3:$WORK/parity/py2/model_germline_indel.json:$WORK/parity/py3/model_germline_indel.json"

        step parity_cp2b 'ptag, max' bash -c '
            set -eu
            R="$1"; P2B="$2"; RT="$3"; T3="$4"; W="$5"; P2LIB="$6"; P2PKGS="$7"; P3PKGS="$8"
            T="$R/scratch/parity/tools"
            OUT="$W/cp2b"; mkdir -p "$OUT"
            for m in "somatic.snv somatic model_somatic_snv somatic_snv_features.csv" \
                     "germline.snv germline model_germline_snv germline_snv_features.csv" \
                     "germline.indel germline model_germline_indel germline_indel_features.csv"; do
                set -- $m
                FSET=$1; SIDE=$2; MODEL=$3; CSV=$4
                env LD_LIBRARY_PATH="$P2LIB" PYTHONPATH="$P2PKGS" \
                    "$P2B" "$T/score_check.py" \
                    --libdir "$RT/$SIDE/lib" \
                    --pkl "$W/parity/py2/$MODEL.pkl" \
                    --csv "$W/parity/py2/$CSV" \
                    --features "$FSET" > "$OUT/${FSET//./_}_py2.tsv"
                env PYTHONPATH="$P3PKGS" python3 "$T/score_check.py" \
                    --libdir "$T3/$SIDE/lib" \
                    --pkl "$W/parity/py3/$MODEL.pkl" \
                    --csv "$W/parity/py3/$CSV" \
                    --features "$FSET" > "$OUT/${FSET//./_}_py3.tsv"
            done
            python3 - "$OUT" <<PY
import sys, os
out = sys.argv[1]
parts = []
for fset in ("somatic_snv", "germline_snv", "germline_indel"):
    def read(p):
        rows = {}
        with open(p) as fid:
            for line in fid:
                line = line.rstrip("\n")
                pos, ref, alt, ptag, ptp = line.split("\t")
                rows[pos] = (ref, alt, ptag, float(ptp))
        return rows
    a = read(os.path.join(out, fset + "_py2.tsv"))
    b = read(os.path.join(out, fset + "_py3.tsv"))
    common = set(a) & set(b)
    ptag_bad = sum(1 for p in common if a[p][2] != b[p][2])
    max_d = max((abs(a[p][3] - b[p][3]) for p in common), default=0.0)
    parts.append("%s: %d/%d ptag, max|dP(TP)| %.4f"
                 % (fset, len(common) - ptag_bad, len(common), max_d))
print(" ; ".join(parts))
PY' _ "$ROOT" "$PY2_BIN" "$REF_TREE" "$PY3_TREE" "$WORK" "$PY2_LIB" "$PY2_PKGS" "$PY3_PKGS"
    else
        step_na parity_cp1 "parity chains not both completed"
        step_na parity_cp2 "parity chains not both completed"
        step_na parity_cp2b "parity chains not both completed"
    fi
fi

# ------------------------------------------------------------ cp3 (C++ runs) layer
if [ $CP3 -eq 1 ]; then
    SNV_M2="$WORK/parity/py2/model_somatic_snv.json"
    SNV_M3="$WORK/parity/py3/model_somatic_snv.json"
    GSNV_M2="$WORK/parity/py2/model_germline_snv.json"
    GSNV_M3="$WORK/parity/py3/model_germline_snv.json"
    GIND_M2="$WORK/parity/py2/model_germline_indel.json"
    GIND_M3="$WORK/parity/py3/model_germline_indel.json"

    cp3_somatic_run() {
        local outdir=$1 model=$2
        mkdir -p "$outdir"
        ( cd "$outdir" && \
          "$INSTALL_DIR/bin/configureStrelkaSomaticWorkflow.py" \
            --tumorBam="$DEMO_DATA/NA12891_demo20.bam" \
            --normalBam="$DEMO_DATA/NA12892_demo20.bam" \
            --referenceFasta="$DEMO_DATA/demo20.fa" \
            --callMemMb=1024 --exome \
            --snvScoringModelFile="$model" \
            --runDir=run && \
          ./run/runWorkflow.py -m local -j 1 -g 4 )
    }
    cp3_germline_run() {
        # genome mode: --exome would force isEVS=False (pre-existing,
        # strelkaGermlineWorkflow.py) and the custom models would not apply
        local outdir=$1 snv=$2 indel=$3
        mkdir -p "$outdir"
        ( cd "$outdir" && \
          "$INSTALL_DIR/bin/configureStrelkaGermlineWorkflow.py" \
            --bam="$DEMO_DATA/NA12891_demo20.bam" \
            --bam="$DEMO_DATA/NA12892_demo20.bam" \
            --referenceFasta="$DEMO_DATA/demo20.fa" \
            --callMemMb=1024 --disableSequenceErrorEstimation \
            --snvScoringModelFile="$snv" \
            --indelScoringModelFile="$indel" \
            --runDir=run && \
          ./run/runWorkflow.py -m local -j 1 -g 4 )
    }

    if has_install && [ -f "$SNV_M2" ] && [ -f "$SNV_M3" ] \
        && [ -f "$GSNV_M2" ] && [ -f "$GIND_M3" ]; then
        step cp3_somatic_m2 - cp3_somatic_run "$WORK/cp3/somatic_m2" "$SNV_M2"
        step cp3_somatic_m3 - cp3_somatic_run "$WORK/cp3/somatic_m3" "$SNV_M3"
        step cp3_germline_m2 - cp3_germline_run "$WORK/cp3/germline_m2" "$GSNV_M2" "$GIND_M2"
        step cp3_germline_m3 - cp3_germline_run "$WORK/cp3/germline_m3" "$GSNV_M3" "$GIND_M3"
        step cp3_check 'identical;' python3 "$QA_TOOLS/cp3_check.py" \
            --somatic-snvs \
                "$WORK/cp3/somatic_m2/run/results/variants/somatic.snvs.vcf.gz" \
                "$WORK/cp3/somatic_m3/run/results/variants/somatic.snvs.vcf.gz" \
            --somatic-indels \
                "$WORK/cp3/somatic_m2/run/results/variants/somatic.indels.vcf.gz" \
                "$WORK/cp3/somatic_m3/run/results/variants/somatic.indels.vcf.gz" \
            --germline \
                "$WORK/cp3/germline_m2/run/results/variants/variants.vcf.gz" \
                "$WORK/cp3/germline_m3/run/results/variants/variants.vcf.gz"
    else
        step_na cp3_somatic_m2 "cp3 prerequisites absent (install tree or parity model JSONs)"
        step_na cp3_somatic_m3 "cp3 prerequisites absent (install tree or parity model JSONs)"
        step_na cp3_germline_m2 "cp3 prerequisites absent (install tree or parity model JSONs)"
        step_na cp3_germline_m3 "cp3 prerequisites absent (install tree or parity model JSONs)"
        step_na cp3_check "cp3 prerequisites absent (install tree or parity model JSONs)"
    fi
fi

# ------------------------------------------------------------- report + env

BOOST=""
BOOST_HDR=$(find "$ROOT/build/bootstrap/boost" -maxdepth 4 -path "*/boost/version.hpp" 2>/dev/null | head -1)
if [ -n "$BOOST_HDR" ]; then
    BOOST=$(grep -m1 "define BOOST_VERSION " "$BOOST_HDR" | awk '{print $3}')
fi
[ -z "$BOOST" ] && BOOST="not materialized (bootstrap at configure time)"

CMAKE_VER=""
if [ -x "$ROOT/build/bootstrap/cmake/bin/cmake" ]; then
    CMAKE_VER=$("$ROOT/build/bootstrap/cmake/bin/cmake" --version 2>/dev/null | head -1)
else
    CMAKE_VER=$(cmake --version 2>/dev/null | head -1 || echo "absent")
fi
GCC_VER=$(gcc --version 2>/dev/null | head -1 || echo "absent")
PY3_VER=$(python3 --version 2>&1)
GIT_BRANCH=$(git -C "$ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")
GIT_DIRTY=$(git -C "$ROOT" status --porcelain 2>/dev/null | wc -l)
PY2_STATE="absent"
has_py2 && PY2_STATE=$(env LD_LIBRARY_PATH="$PY2_LIB" "$PY2_BIN" --version 2>&1)
PY3ENV_STATE="absent"
has_py3env && PY3ENV_STATE="present"
INSTALL_STATE="absent"
has_install && INSTALL_STATE="present"

python3 "$QA_TOOLS/qa_report.py" \
    --out "$QA_DIR/evidence" \
    --stamp "$STAMP" \
    --level "$LEVEL" \
    --git-sha "$GITSHA" \
    --steps "$STEPS" \
    --env "git_sha=$GITSHA" \
    --env "git_branch=$GIT_BRANCH" \
    --env "git_dirty_files=$GIT_DIRTY" \
    --env "python3=$PY3_VER" \
    --env "gcc=$GCC_VER" \
    --env "cmake=$CMAKE_VER" \
    --env "boost=$BOOST" \
    --env "install_tree=$INSTALL_STATE" \
    --env "toolchain_py2=$PY2_STATE" \
    --env "toolchain_py3_env=$PY3ENV_STATE" \
    --env "py2_ref_tree_sha=$PY2_TREE_SHA" \
    --input "scratch/parity/checksums.sha256=$ROOT/scratch/parity/checksums.sha256" \
    --input "requirements-training.txt=$PY3_TREE/requirements-training.txt" \
    --input "src/demo/expectedResults=$ROOT/src/demo/expectedResults" \
    --input "src/demo/expectedResultsGermline=$ROOT/src/demo/expectedResultsGermline" \
    --input "redist=$ROOT/redist"
rc=$?

echo
echo "qa gate ($LEVEL): report in $QA_DIR/evidence/ (rc=$rc: 0=PASS 1=FAIL)"
exit $rc
