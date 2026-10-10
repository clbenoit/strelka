# Strelka modernization merge log — 2026-10-10

Integration of the two OpenSpec modernization changes into the `python3` branch,
sequentially, with `--no-ff` merge commits so the branch topology stays visible
in the repository history.

## Changes integrated

| # | OpenSpec change | Branch | Commits | Scope |
|---|---|---|---|---|
| 1 | `strelka-py3-runtime` (archived 2026-10-09) | `openspec/strelka-py3-runtime` | `84ebde6f`, `9b96ab4f` | Python 2 → 3.8+ runtime (workflow bin/lib/libexec, vendored pyflow, CMake `PythonInterp 3`) |
| 2 | `strelka-redist-modern` (archived 2026-10-09) | `openspec/strelka-redist-modern` | `4a6983c2` | Embedded toolchain: htslib/samtools 1.24, Boost 1.92.0 subset, CMake 3.31.12, gcc 14 fixes |

Both branches cut from `python3` @ `4f2a99a4` (2026-10-09), developed in isolated
worktrees, validated independently, then merged sequentially:
**1) py3-runtime → 2) redist-modern** (the redist E2E validation ran with the py3
runtime overlaid — the overlay becomes permanent with this merge).

## Resulting topology

```
*   5a77770f Merge openspec/strelka-redist-modern (M2)
|\
| * 4a6983c2 redist: htslib/samtools 1.24, Boost 1.92.0, CMake 3.31.12
* | 35fa60d5 Merge openspec/strelka-py3-runtime (M1)
|\ \
| |/
|/|
| * 9b96ab4f py3: fix 3 runtime bugs caught by E2E demos
| * 84ebde6f py3: modernize runtime to Python 3.8+
|/
* 4f2a99a4 (base, python3 = v2.9.x)
```

- M1 = `35fa60d5` (fast-forwardable at the base; `--no-ff` keeps the branch line)
- M2 = `5a77770f` (necessarily a true merge — M1 already moved `python3`)

## Conflict resolution record

**Textual conflicts: 0** — verified with `git merge-tree --write-tree` in both
orderings before merging. Three files were touched by both branches and merged
automatically; the merged content of each was manually reviewed (union correct):

| Shared file | py3-runtime side | redist-modern side | Resolution |
|---|---|---|---|
| `CMakeLists.txt` | `find_package(PythonInterp 3)`, min python 3.8 | `cmake_minimum_required(3.10)` | auto-merge, disjoint hunks — both sets present |
| `docs/userGuide/installation.md` | python 2.6+ → 3.8+ | CMake 3.10+ / Boost 1.92.0+ prerequisites | auto-merge, disjoint sections |
| `docs/developerGuide/README.md` | "python 2 (2.6+)" → "python 3 (3.8+)" on L194 | Boost 1.58 → 1.92.0 notes (L107/111) | auto-merge; **one semantic gap below** |

**Semantic conflict: 1** — after M2, the merged `docs/developerGuide/README.md` L194
still read "Supported project languages are **C++11** …": each branch alone was
correct (py3 branch kept the C++11 toolchain, redist branch kept the python2
wording), but the union predates the C++14 bump that redist-modern required
(Boost.Math 1.92 headers need C++14). Resolved by the post-merge follow-up
commit: "C++11" → "C++14" on that line (the VS2013 note at L163 was left as-is:
Windows/VS specifics, out of scope of the container-validated build).

## Build saga (redist-modern, for the record)

The gcc 14 / toolchain build converged in 8 make iterations, each fix minimal and
kept at the break point (no refactoring):

1. `bjam` → `b2` (boost 1.92 bootstrap no longer creates `bjam`)
2. cmake version check → prefix match (multi-line output since cmake 3.0)
3. `tools/boost_install` kept in the subset (b2 install rules live there since ~1.77)
4. subset closure re-computed on all jam files (build.jam + legacy Jamfile.v2 stubs +
   implicit `/boost//headers` dependency) → 64 dirs incl. `headers`, `align`, `scope`,
   `io`, `phoenix`, `numeric/conversion`; validated with `b2 -n --debug-loading`
5. C++ standard 11 → 14 (`cxxConfigure.cmake`) — Boost.Math 1.92 requirement
6. rapidjson 1.1.0 one-line re-pack (`GenericStringRef::operator=` ill-formed in C++14)
7. htslib `--without-libdeflate` (system libdeflate auto-detected but never linked)
8. misc: `stringer<T>()` template-id ctor, `[[noreturn]]` on qscore_cache error fns,
   UB read padded in `IndelModelProduction.cpp` (theta-locked minimizer workspace),
   `-Wno-maybe-uninitialized` scoped to the two `*_sim_test.cpp` (boost random
   gcc 14 false positive), missing `#include <limits>` in `IndelErrorModel.hh`

The E2E workflow demos (first full run of the py3 runtime path) additionally caught
3 conversion bugs in the py3-runtime branch, fixed in its follow-up commit `9b96ab4f`:

- `configureUtil.py`: spurious `.items()` on `config.items(section)`
- `checkChromSet.py`: 2× `subprocess.Popen` needed `universal_newlines` (bytes vs str)
- `workflowUtil.py`: py2 int divisions left as floats → `range(float)` TypeError

## Validation on the merged tree

Run in an isolated worktree on the merged head (fresh configure, bootstrapped
cmake + boost stages reused):

| Check | Result |
|---|---|
| `configure` + `make -j16` | exit 0, first try (no fix needed on the union) |
| `ctest` | **9/9** — the 9th suite (`strelka_unit_test_python_lib`) is new: python unit tests now run under ctest because `PythonInterp 3` is found |
| `make install` | exit 0 (`configBuildTimeInfo.py` now generated at build stage) |
| `runStrelkaGermlineWorkflowDemo.bash` | exit 0 — "Demo/verification successfully completed" (native python3, no overlay) |
| `runStrelkaSomaticWorkflowDemo.bash` | exit 0 — "No differences between expected and computed results" |
| Effective versions | samtools 1.24, htslib 1.24, Boost 1.92.0 (109200), cmake 3.31.12 |

## Post-merge state

- `python3` pushed to origin (this document is part of the follow-up commit)
- OpenSpec store: both changes archived, specs synced (`specs/strelka-python3-runtime`,
  `specs/strelka-redist-modern`); post-archive integration note recorded in
  `archive/2026-10-09-strelka-redist-modern/tasks.md`
- Branch cleanup pending user validation: `git worktree remove` +
  `git branch -d openspec/strelka-py3-runtime openspec/strelka-redist-modern`
