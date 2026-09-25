# Repository identity and branch checks

This checkout is for **skilledwolf/ALPS**, even when its directory is named
`ALPS_upstream`. The directory name does not identify the Git remote or branch.

- `origin` is `https://github.com/skilledwolf/ALPS.git`, the primary fetch,
  tracking, and push remote.
- `alpsim` is `https://github.com/ALPSim/ALPS.git`, the original project. Use it
  as a task's target only when the user explicitly requests that project.
- The CMake/SDK/CI overhaul currently targets
  `codex/windows-ci-runner-setup`, tracking
  `origin/codex/windows-ci-runner-setup`. An explicit user instruction to use
  another branch takes precedence.
- The sibling `../ALPS` checkout is a separate repository (`skilledwolf/alpsx`);
  do not substitute it for this checkout.

## Before reviewing or editing code

Run these checks at the start of each task and after switching checkouts:

```sh
git rev-parse --show-toplevel
git remote -v
git status --short --branch
git branch --show-current
git rev-parse --abbrev-ref '@{upstream}'
git log -1 --format='%h %s'
git rev-list --left-right --count 'HEAD...@{upstream}'
```

Establish the intended repository, branch, and commit from the user's request
and these results before reading code for an assessment or making changes.
State the repository and branch in the initial progress update. Confirmation
of a directory alone is not confirmation of its currently checked-out branch.

When the user requests the latest feature branch state, fetch `origin` and
compare against the intended remote branch. A stale local tracking ref or the
most recently timestamped branch is not sufficient evidence of the target.

If the checkout does not match the intended branch, resolve that mismatch
before doing the work. Preserve unrelated changes and use the user's existing
branch instructions; ask for clarification only when the target is genuinely
ambiguous or changes conflict. Do not review an old release branch and present
its findings or test results as applying to the current feature branch.

## Build overhaul and changes

The maintainer explicitly requires no legacy build compatibility for this
overhaul. Remove obsolete build interfaces instead of preserving aliases or
fallbacks solely for compatibility. Keep the intended scientific functionality
and validate changes against the actual target branch.

Keep unrelated work and untracked review artifacts intact. Never reset or
clean the whole checkout to switch branches. Do not push changes unless the
user explicitly requests a push.
