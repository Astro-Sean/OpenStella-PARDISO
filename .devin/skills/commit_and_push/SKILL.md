---
name: commit_and_push
description: Commit each changed file in the OpenStella_parallel tree individually with a short humanized message, then push to Astro-Sean/OpenStella-PARDISO
allowed-tools:
  - exec
  - read
  - grep
  - glob
  - skill
---

Commit and push the OpenStella_parallel working tree. One commit per file, short messages, then push to `https://github.com/Astro-Sean/OpenStella-PARDISO`.

## Rules

- One commit per modified file. Stage a single file, commit it, move to the next. Exception: a source change and the documentation note that explains it may share a commit when the note only makes sense with that change.
- Commit messages are one line, imperative mood, under 72 characters. State what changed in that file. Prefix with the routine or area when it helps (`strad: ...`, `eve: ...`, `build: ...`), matching `git log`. No body text unless a change is genuinely hard to summarize.
- Write every message in plain human prose, then run it through the humanizer skill (invoke `humanizer`, or apply its rules if it is not available): no not-X-but-Y contrasts, no dashes as connectors, no inflated words, no staging, no stock AI words. Plain ASCII only.
- Never mention AI tools or assistance. No `Co-Authored-By` trailers, no "Generated with" lines, no tool names. Commits are authored by the maintainer only. Do not change git config.
- Never commit build artifacts, backups, or run products. Do not `git add -f`. Leave these files uncommitted and report them: `.built_*`, `*.o`, `*.mod`, `*.lst`, `*.out`, `*.err`, `*.log`, `*.exe`, `*.exe.*`, `*.bak`, `*.bak_*`, `*.backup`, `*~`, `*.swp`, `core.*`, `strad.1`, `obj/*.inc`, `bin/` executables, `run/` executables, `pardiso_dump_*.dat`, `pardiso_badmat_*.dat`, and files ignored by `.gitignore`.
- Only commit source changes that build clean. If a `.f`, `.f90`, `.inc`, `Makefile`, `*.mak`, or `build*.sh` change has not been through `./build_parallel.sh` with zero errors since it was edited, stop and say so. Comment-only and doc-only edits do not need a rebuild.
- Push only to the OpenStella-PARDISO remote. `origin` must be `https://github.com/Astro-Sean/OpenStella-PARDISO`, `https://github.com/Astro-Sean/OpenStella-PARDISO.git`, or `git@github.com:Astro-Sean/OpenStella-PARDISO.git`. If `origin` points anywhere else, stop and ask before changing remotes.
- Push the current branch to its upstream branch on `origin`. If the branch has no upstream, or the upstream is not on `origin`, stop and ask.
- Do not rewrite history, amend published commits, or force-push.
- If nothing is modified, say so and stop. Do not commit.

## Workflow

1. Inspect the tree and confirm the push target:

   ```bash
   git status --short
   git diff
   git log --oneline -8
   git remote get-url origin
   git branch --show-current
   git rev-parse --abbrev-ref @{upstream}
   ```

   The remote must be OpenStella-PARDISO and the upstream must be on `origin`. If either is wrong, stop and ask. Read the `git log` output to match the existing message style.

2. For each modified or new source file, read its diff, stage it alone, and commit:

   ```bash
   git add <file>
   git commit -m "<short imperative message>"
   ```

   Repeat until only ignored or excluded files remain.

3. Push to the branch's upstream on `origin`:

   ```bash
   git push
   ```

   If the push is rejected as non-fast-forward, stop and report. Do not force.

4. Report the commit hashes and files in a short list, plus any files left uncommitted and why.
