#!/usr/bin/env python3
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


"""Non-regression evidence report writer (change strelka-qa-noregression).

Consumes the gate's manifest (steps TSV + environment facts + pinned inputs)
and writes ONE append-only, self-hashed evidence report:

    <out>/<stamp>.json   machine-readable record (sha256 of everything)
    <out>/<stamp>.md     human-readable summary, carries the .json sha256

Rules (qa/README.md):
  - append-only: existing report files are never modified; a stamp collision
    gets a -1, -2, ... suffix
  - self-hashed: the .md carries the sha256 of the .json, so a record is
    verifiable on its own
  - explicit N/A: rc -1 = step explicitly skipped (e.g. toolchain absent),
    never silent
  - self-contained: each step carries a one-line `result` (key numbers) and,
    on failure, a tail of its log — the report does not depend on temp files

Steps TSV columns (tab-separated, one line per step):
    name  rc  logpath  result
"""

import argparse
import datetime
import hashlib
import json
import os
import sys


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fid:
        for chunk in iter(lambda: fid.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_tree(path):
    """Aggregate hash of a directory: per-file sha256, manifest sorted by
    repo-relative path, sha256 of the manifest (stable walk order)."""
    h = hashlib.sha256()
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for name in sorted(files):
            full = os.path.join(root, name)
            rel = os.path.relpath(full, path)
            h.update(rel.encode("utf-8"))
            h.update(sha256_file(full).encode("ascii"))
    return h.hexdigest()


def tail(path, n=15):
    with open(path) as fid:
        lines = fid.read().splitlines()
    return "\n".join(lines[-n:])


def parse_steps(path):
    steps = []
    with open(path) as fid:
        for line in fid:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            while len(parts) < 4:
                parts.append("")
            name, rc, log, result = parts[0], int(parts[1]), parts[2], parts[3]
            step = {"name": name, "rc": rc, "result": result}
            step["log_sha256"] = (
                sha256_file(log) if os.path.isfile(log) else None)
            if rc > 0 and os.path.isfile(log):
                step["log_tail"] = tail(log)
            steps.append(step)
    return steps


def free_name(out_dir, stamp):
    base = os.path.join(out_dir, stamp)
    for i in range(10):
        suffix = "" if i == 0 else "-%d" % i
        if not (os.path.exists(base + suffix + ".json")
                or os.path.exists(base + suffix + ".md")):
            return base + suffix
    raise Exception("stamp collision unresolved: %s" % stamp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="evidence directory")
    parser.add_argument("--stamp", required=True, help="<UTC>-<gitsha7>")
    parser.add_argument("--level", required=True,
                        choices=["fast", "full", "cp3"],
                        help="gate level actually executed")
    parser.add_argument("--git-sha", required=True)
    parser.add_argument("--steps", required=True,
                        help="TSV: name<TAB>rc<TAB>logpath<TAB>result")
    parser.add_argument("--env", action="append", default=[],
                        help="environment fact KEY=VALUE (repeatable)")
    parser.add_argument("--input", action="append", default=[],
                        help="pinned input LABEL=<path> (repeatable; "
                             "file = sha256, directory = aggregate hash)")
    args = parser.parse_args()

    steps = parse_steps(args.steps)
    env = {}
    for item in args.env:
        key, _, value = item.partition("=")
        env[key] = value
    inputs = {}
    for item in args.input:
        label, _, path = item.partition("=")
        if os.path.isfile(path):
            inputs[label] = sha256_file(path)
        elif os.path.isdir(path):
            inputs[label] = sha256_tree(path) + " (aggregate)"
        else:
            inputs[label] = "MISSING (%s)" % path
            print("warning: pinned input missing: %s" % path, file=sys.stderr)

    fails = [s for s in steps if s["rc"] > 0]
    record = {
        "stamp": args.stamp,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc)
        .strftime("%Y-%m-%dT%H:%M:%SZ"),
        "level": args.level,
        "git_sha": args.git_sha,
        "result": "FAIL" if fails else "PASS",
        "steps": steps,
        "steps_na": [s["name"] for s in steps if s["rc"] == -1],
        "steps_failed": [s["name"] for s in fails],
        "environment": env,
        "pinned_inputs": inputs,
    }

    os.makedirs(args.out, exist_ok=True)
    base = free_name(args.out, args.stamp)
    json_path = base + ".json"
    with open(json_path, "w") as fid:
        json.dump(record, fid, indent=2, sort_keys=True)
        fid.write("\n")
    json_hash = sha256_file(json_path)

    lines = []
    lines.append("# Strelka non-regression evidence — %s" % args.stamp)
    lines.append("")
    lines.append("- result: **%s** (level: %s)" % (record["result"],
                                                   args.level))
    lines.append("- git: %s" % args.git_sha)
    lines.append("- generated: %s" % record["generated_utc"])
    lines.append("- report sha256 (json): `%s`" % json_hash)
    lines.append("")
    lines.append("## Steps")
    lines.append("")
    lines.append("| step | rc | verdict | result | log sha256 |")
    lines.append("|---|---|---|---|---|")
    for s in steps:
        verdict = ("PASS" if s["rc"] == 0 else
                   "N/A (explicit)" if s["rc"] == -1 else "FAIL")
        log = s["log_sha256"] or "-"
        result = s["result"].replace("|", "\\|") or "-"
        lines.append("| %s | %d | %s | %s | `%s` |"
                     % (s["name"], s["rc"], verdict, result, log))
    if record["steps_na"]:
        lines.append("")
        lines.append("Explicitly skipped (consigné, jamais silencieux): %s"
                     % ", ".join("`%s`" % n for n in record["steps_na"]))
    if fails:
        lines.append("")
        lines.append("## Failures (log tails)")
        lines.append("")
        for s in fails:
            lines.append("### %s (rc=%d)" % (s["name"], s["rc"]))
            lines.append("")
            lines.append("```")
            lines.append(s.get("log_tail", "(log unavailable)"))
            lines.append("```")
            lines.append("")
    lines.append("## Environment (inventaire)")
    lines.append("")
    for key in sorted(env):
        lines.append("- %s: %s" % (key, env[key]))
    lines.append("")
    lines.append("## Pinned inputs (sha256)")
    lines.append("")
    for label in sorted(inputs):
        lines.append("- %s: `%s`" % (label, inputs[label]))
    md_path = base + ".md"
    with open(md_path, "w") as fid:
        fid.write("\n".join(lines) + "\n")

    print("evidence report: %s (%s)" % (md_path, record["result"]))
    if fails:
        sys.exit(1)


if __name__ == "__main__":
    main()
