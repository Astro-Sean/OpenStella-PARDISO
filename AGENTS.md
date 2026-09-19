# AGENTS.md

## Commit conventions

- Do not add AI attribution to commits. No "Generated with" lines and no
  "Co-Authored-By" trailers for Devin or any other tool. Commits are authored
  by the repository owner only.
- Write commit messages that explain why a change was made, not just what
  changed.

## Build

- `bash build_parallel.sh <NZON>` builds `bin/xstella.exe` with Intel Fortran
  and MKL. `NZON` is the radial zone count and sets `Mzon` in the generated
  `src/zone.inc`.
- The solver path uses pre-generated `.f` sources in `strad/` and `src/`
  directly; trefor preprocessing is not needed.
- `src/zone.inc`, `src/opacity.inc`, `bin/`, and `*.o`/`*.mod` are generated
  build products and are gitignored. Do not commit them.
