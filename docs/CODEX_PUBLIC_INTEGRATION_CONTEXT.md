# Codex context: porting public SAMsrcV5 changes

Load this entire file before integrating code that was developed from the
original SAMsrcV5 public codebase. Treat the requirements below as project
instructions for the integration task.

## Reusable invocation

Use this context with a request such as:

```text
Read docs/CODEX_PUBLIC_INTEGRATION_CONTEXT.md completely and follow it.
Incoming source ref: <branch, tag, or commit>
Requested scope: <optional subset; otherwise all changes from the baseline>

Port the incoming behavior into the currently checked-out portable branch.
Do not merge or cherry-pick the legacy implementation directly.
```

If the incoming ref or requested scope is missing and cannot be derived from
the task, ask for it before changing code. Resolve every ref to a full commit
SHA and report that SHA in the final result.

## Fixed history and objective

- The immutable original-code baseline is
  `fbf259de66adfa2dd6950839c986c145d68e58ee` (`fbf259d`, the public branch at
  the start of the portable conversion).
- Never replace or advance this baseline to a newer public commit.
- The portable target is the currently checked-out branch containing
  `pyproject.toml`, `CMakeLists.txt`, `src/samsrcv5`, and the split SAM
  input/output implementation. At the time this context was written, that
  branch is `convert2multiOS`.
- The incoming change set is always defined first as the semantic difference
  between the fixed baseline and the supplied incoming source ref. Do not use
  a diff against the portable branch as a substitute: that would mix incoming
  functionality with the packaging and portability conversion.
- The objective is to preserve the incoming behavior while expressing it in
  the portable branch's current architecture: pip-installable platform wheels,
  bundled native programs and required libraries, package-relative resources,
  cross-platform behavior, and independent `-i_SAMdir`/`-o_SAMdir` roots.

## Non-negotiable working rules

1. Inspect before editing. Preserve unrelated user changes and do not overwrite
   or reformat files outside the incoming change's scope.
2. Do not merge, rebase, cherry-pick, or blindly apply the legacy diff to the
   portable target. Reimplement each relevant change against the current code.
3. Do not assume every baseline-to-incoming hunk still needs work. Public
   changes may have been ported previously or superseded in the portable
   branch. Prove the disposition of every relevant hunk.
4. Preserve both build paths where they remain supported: the pip/CMake wheel
   build is the distribution path, while the legacy Make build remains a Unix
   developer path.
5. Do not weaken, remove, or silently bypass portability, packaging, licensing,
   CLI, or test behavior just to make an incoming patch apply.
6. Do not commit, push, fetch, switch branches, or alter Git refs unless the
   user's request authorizes that action. Read-only Git inspection is expected.
7. Keep the implementation focused on the incoming behavior. Do not fold in
   unrelated cleanup unless it is necessary for a correct portable port.

## Phase 1: establish the change set

Start with read-only repository checks. Use task-specific variable names rather
than generic shell variables:

```sh
SAM_LEGACY_BASE=fbf259de66adfa2dd6950839c986c145d68e58ee
SAM_INCOMING_REF=<incoming-ref>

git status --short
git branch --show-current
git cat-file -e "$SAM_LEGACY_BASE^{commit}"
git rev-parse --verify "$SAM_INCOMING_REF^{commit}"
git merge-base --is-ancestor "$SAM_LEGACY_BASE" "$SAM_INCOMING_REF"
git log --reverse --format='%H%x09%s' "$SAM_LEGACY_BASE..$SAM_INCOMING_REF"
git diff --find-renames --find-copies --stat \
  "$SAM_LEGACY_BASE" "$SAM_INCOMING_REF"
git diff --find-renames --find-copies --name-status \
  "$SAM_LEGACY_BASE" "$SAM_INCOMING_REF"
git diff --find-renames --find-copies \
  "$SAM_LEGACY_BASE" "$SAM_INCOMING_REF" -- <relevant paths>
```

The ancestry check must succeed. If it does not, stop and explain that the
incoming ref is not based on the required original commit; ask the user to
identify the correct source ref or explicitly redefine the requested scope.

Read the current portable versions of every affected file and their nearby
tests before deciding how to port a hunk. Also inspect the portable conversion
history when it explains why the target differs:

```sh
git log --oneline --decorate \
  fbf259de66adfa2dd6950839c986c145d68e58ee..HEAD
git diff --find-renames \
  fbf259de66adfa2dd6950839c986c145d68e58ee HEAD -- <affected paths>
```

For each incoming file or coherent behavior, maintain a working inventory with
one of these dispositions:

- **Port:** behavior is absent and must be implemented in the portable style.
- **Already present:** equivalent behavior exists; cite the target code/test.
- **Superseded:** the portable implementation intentionally solves the same
  need differently; explain why no legacy hunk should be copied.
- **Not applicable:** generated, build-only, obsolete, or outside an explicitly
  requested scope; give a concrete reason.
- **Blocked:** intent or required dependency cannot be determined safely; ask
  the user rather than guessing.

For cumulative public refs, the fixed-baseline diff may include changes ported
in an earlier session. Always perform this reconciliation so those changes are
not duplicated or reverted.

## Phase 2: transform changes into the portable architecture

### Native libraries and shared C code

- Add new reusable C sources to the appropriate legacy `Libs/*` Make build and
  to the explicit `SAM_LIBRARY_SOURCES` list in `CMakeLists.txt` so they become
  part of `samcore` and the installed native programs.
- Add or update public declarations in `include/` and keep callers consistent.
- Retain C17 compatibility. Avoid compiler extensions unless they are already
  deliberately supported on every wheel platform.
- Account for Linux, macOS, and native Windows/MinGW-w64. Reuse
  `include/sam_platform.h` for narrowly scoped compatibility shims rather than
  scattering incompatible platform assumptions through scientific code.
- Open binary scientific formats with binary modes (`rb`, `wb`, or equivalent),
  which is required on Windows.
- Do not assume POSIX-only path syntax, `/tmp`, `PWD`, `HOME`, Unix executable
  suffixes, or availability of functions absent from MinGW. Use existing
  portable helpers and the `USERPROFILE`, `TEMP`, and executable-suffix
  conventions where applicable.
- Preserve the numerical algorithm and file format unless the incoming change
  explicitly changes them. Add regression evidence for any numerical change.

### Native executables

A native command is not integrated merely because its source compiles with
Make. For every new or newly supported executable:

1. Add an explicit CMake target linked to `samcore` and any required libraries.
2. Install the executable into `samsrcv5/_bin` through CMake.
3. Add a launcher function in `src/samsrcv5/launcher.py` using the existing
   `_native` mechanism, including Windows `.exe` and repaired-DLL lookup.
4. Add the public command name to `[project.scripts]` in `pyproject.toml`.
5. Preserve established command names and intentional `.py` aliases.
6. Add source-build CLI tests and installed-wheel tests for command presence,
   help/error behavior, and at least one meaningful behavior where feasible.

Do not expose a command only through a repository-local `bin/` directory.

### Python code and legacy scripts

- Prefer a maintained module under `src/samsrcv5` for new Python functionality.
  Keep a script under packaged `_legacy` only when its existing script-style
  imports and behavior make that the safer compatibility choice.
- Register user-facing commands in `[project.scripts]`; do not depend on Make
  substituting `@@libdir@@` for wheel installations.
- Resolve bundled files with `importlib.resources` or paths based on the
  installed package, never the repository checkout or current working
  directory.
- Use `pathlib`, `os.pathsep`, `sysconfig`, and platform-aware subprocess
  handling. Do not hard-code `/`, `:` as a path-list separator, or Unix command
  locations.
- Declare runtime Python dependencies in `pyproject.toml` with the narrowest
  justified compatibility constraints. Do not add an undeclared import.
- Keep optional host tools such as AFNI, FreeSurfer, qhull, GTK/PyGObject, and
  Tk as runtime checks when they are not link-time wheel dependencies. Help
  output and unrelated commands should continue to work when optional tools
  are absent.

### Package data

- Explicitly install new templates, UI files, atlases, or other runtime assets
  into `samsrcv5/data` in `CMakeLists.txt`.
- Access those assets through the installed package. A wheel test must assert
  that each required asset exists and, when practical, can be opened.
- Include all required source paths in the sdist configuration. Verify both a
  source archive and a wheel rather than relying on checkout-only tests.

### Native and Python dependencies

- Wheels must install without requiring the user to provide a compiler, GSL,
  FFTW, or other intended bundled native runtime libraries.
- For a new native dependency, decide explicitly whether it is statically
  linked, repaired into the wheel, or remains an optional external workflow
  tool. Do not leave this decision to the installer.
- Update dependency download/build scripts, hashes, CI repair/inspection, source
  distribution contents, license files, `THIRD_PARTY_NOTICES.md`, and
  corresponding-source artifacts when required by the dependency's license.
- Ensure dependency handling works for Linux x86_64/aarch64, macOS x86_64/arm64,
  and Windows AMD64. A successful local Linux link is not sufficient evidence.

## SAM input and output directory contract

Every incoming access to a SAM product must be classified as a read or a write.
Transform hard-coded forms such as `%d/SAM`, `<DSpath>/SAM`, and paths derived
directly from the dataset directory according to this contract.

### Required semantics

- `-i_SAMdir <SAMDIR>` and parameter-file/long-option
  `InputSAMDirectory <SAMDIR>` select the root for existing SAM products.
- `-o_SAMdir <SAMDIR>` and parameter-file/long-option
  `OutputSAMDirectory <SAMDIR>` select the root for newly created SAM products,
  run parameters, and logs.
- Each side defaults independently to `<dataset>/SAM` when omitted. Supplying
  only one option must not redirect the other side.
- An option names the SAM root itself; do not append another `SAM` component.
- Validate required input roots and files without creating them.
- Create output roots and missing parents with the existing recursive portable
  helper before writing.
- Explicit `ImageDirectory` behavior remains an override for final images where
  the current command supports it. It must not silently redirect other output
  products or logs.
- Preserve parameter-file behavior and command-line precedence. Parameter files
  remain optional for the commands that currently support command-line-only
  operation, including the `cmdline` naming fallback.

### Parser choice

- Programs using the shared SAM parameter system must register the standard
  parameters, populate `PARMINFO.InputSAMDirectory` and
  `PARMINFO.OutputSAMDirectory`, and resolve roots with
  `GetSAMPath(..., SAM_INPUT)` or `GetSAMPath(..., SAM_OUTPUT)`.
- Older standalone `getopt` programs must call `parse_samdir_args` before
  `getopt`/`getopt_long` sees the multi-character single-dash flags, document
  both flags in usage, apply independent defaults, and use portable input
  validation/output creation.
- Do not introduce a third implementation of the same parsing or fallback
  logic when an existing shared helper can be extended safely.

### Read/write audit examples

- Covariance, weights, transforms, targets, masks, and other pre-existing SAM
  products are read from the input root.
- New covariance directories, weights, forward solutions, images, noise files,
  parameter reconstructions, and logs are written under the output root unless
  a documented explicit output option overrides that particular product.
- A command that both consumes and produces a same-named product directory must
  still use separate roots. Never infer that equal product names imply equal
  directories.

Add tests that run with distinct input and output roots and assert that no
dataset-local `SAM` directory is created. Also retain a test proving that both
omitted options preserve the legacy `<dataset>/SAM` default.

## Documentation and compatibility

- Update user-facing installation, command, and option documentation when the
  incoming behavior changes them.
- Preserve documented public command names and file formats unless the user has
  explicitly authorized a breaking change.
- Keep `docs/PORTING.md`, README installation guidance, and the actual CMake,
  `pyproject.toml`, launchers, and CI configuration consistent.
- Treat Makefiles from the incoming ref as evidence of sources, generated data,
  dependencies, and install behavior. Translate that evidence into CMake and
  wheel metadata; do not copy Make-only installation assumptions unchanged.

## Verification ladder

Choose checks in proportion to the change, but do not omit a relevant layer.
Report commands that were not run and why.

1. **Static checks:** inspect the final diff, run `git diff --check`, search for
   stale hard-coded SAM paths and source-tree resource assumptions, and verify
   all new files are represented in CMake/sdist/package metadata.
2. **Targeted tests:** run or add focused C/Python tests for every incoming bug
   fix or new behavior, including failure paths.
3. **Host unit and CLI tests:** run `make test-unit` when native source/build
   prerequisites are available.
4. **CMake/package build:** configure and build the CMake project, then build a
   wheel with the repository's documented dependency prefix. Do not treat the
   legacy Make build alone as package validation.
5. **Installed-wheel tests:** install the built wheel into a clean isolated
   environment and run `test/wheel` from outside the checkout so local sources
   cannot shadow the installed package.
6. **AFNI/CTF integration:** run `make test-integration` for changes affecting
   SAM products, directory routing, scientific I/O, native algorithms, AFNI,
   CTF data, orthohull, or pipeline commands. Use `make test-slow` only when the
   full-volume behavior is relevant or specifically requested.
7. **Cross-platform wheels:** changes to C/CMake, dependencies, launchers,
   filesystem behavior, binary I/O, or package contents require the GitHub
   wheel matrix to validate Linux, macOS, and Windows. Note this explicitly if
   it cannot be reproduced locally.

Before any Python package installation or mutation, inspect the intended
interpreter with read-only checks such as:

```sh
command -v python3
python3 -c 'import sys; print(sys.executable); print(sys.prefix)'
printf '%s\n' "${CONDA_PREFIX-}"
```

Never install, update, or remove packages in a Conda-managed environment
without explicit user approval for the exact environment, command, and
packages. Prefer a project-local `.venv` created with a suitable non-Conda
Python, invoke its Python and pip explicitly, and confirm `.venv` remains
ignored by Git. Do not silently fall back to a user or system installation.

## Definition of done

The integration is complete only when all of the following are true:

- The exact fixed-baseline-to-incoming diff has been inventoried.
- Every relevant incoming behavior has a documented disposition.
- Ported behavior follows the current portable architecture rather than
  restoring original installation or path assumptions.
- New C library code and native executables are included in the wheel build.
- New Python commands and data work from an installed wheel outside the source
  tree.
- All SAM reads and writes honor independent input/output roots and legacy
  defaults.
- New dependencies and bundled code have complete platform, source, and license
  handling.
- Relevant tests pass, or unavailable test layers and their remaining risk are
  clearly reported.
- `git diff --check` passes and unrelated user changes remain untouched.

## Required final report

Conclude each integration with a concise report containing:

1. The fixed baseline SHA, exact incoming SHA, and target branch/HEAD used.
2. A summary of incoming behaviors and their dispositions: ported, already
   present, superseded, not applicable, or blocked.
3. The portable transformations made, especially packaging, dependencies,
   resources, platform handling, and SAM input/output routing.
4. Tests and build checks run with outcomes.
5. Checks not run, why they were unavailable, and which CI jobs must supply the
   remaining evidence.
6. Any compatibility, licensing, numerical, or follow-up concerns.

Do not describe a legacy patch as integrated merely because it applies or
compiles. The measure of completion is equivalent intended behavior delivered
through the pip-installable, cross-platform SAMsrcV5 package.
