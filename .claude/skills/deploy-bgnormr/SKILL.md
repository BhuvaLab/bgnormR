---
name: deploy-bgnormr
description: Deploy / release / publish bgnormR to GitHub. Runs the full pre-push ritual — document, update vignette, code review, R CMD check, BiocCheck — fixing issues, then bumps the version and commits & pushes. Use when asked to deploy, release, ship, publish, or "run the checks before pushing" bgnormR, or prepare a Bioconductor push.
---

# Deploy bgnormR

bgnormR is a Bioconductor-track package (0.99.z). "Deploying" means running the
full pre-push ritual so GitHub CI
([.github/workflows/check-bioc.yml](.github/workflows/check-bioc.yml)) stays
green and the package stays clean for Bioconductor review. Work through the
steps **in order**; do not skip to the commit.

The mechanical check phases are wrapped by the helper
[.claude/skills/deploy-bgnormr/deploy.R](deploy.R):

```bash
Rscript .claude/skills/deploy-bgnormr/deploy.R <document|check|bioccheck|all>
```

It writes logs to `$DEPLOY_OUT` (default: a tempdir it prints on start) and
exits non-zero when a phase surfaces a **blocking** problem (check errors or
warnings, real BiocCheck errors/warnings). Set `DEPLOY_OUT` to a stable path so
you can read the logs:

```bash
export DEPLOY_OUT=/tmp/bgnorm-deploy
```

All paths below are relative to the repo root (`bgnormR/`).

## Step 0 — Preflight

Confirm the push target with the user. Default to the current branch:

```bash
git rev-parse --abbrev-ref HEAD && git status --short
```

Most runs push the branch already checked out — **confirm it each time** before
the final push. Note what's already modified so you can describe it in the
commit later.

## Step 1 — Document

```bash
Rscript .claude/skills/deploy-bgnormr/deploy.R document
```

Regenerates `man/*.Rd` and `NAMESPACE` from roxygen comments. Review the diff —
`document()` **rewrites `NAMESPACE`**, so a stale `@export`/`@importFrom` shows
up here.

## Step 2 — Update the vignette (judgment step)

Only if a **significant new feature** shipped. Compare exported functions in
[NAMESPACE](NAMESPACE) and recent commits (`git log --oneline -15`) against
[vignettes/bgnormR.Rmd](vignettes/bgnormR.Rmd). If a notable new exported
function or capability isn't demonstrated, add a short section. Skip silently if
nothing significant changed — do not pad the vignette.

If you edit it, confirm it still knits before moving on (a broken vignette fails
`check` in step 4 anyway):

```bash
Rscript -e 'devtools::build_vignettes()'
```

## Step 3 — Code review

Invoke the **`/code-review`** skill on the working diff, triage the findings, and
apply fixes. (You run the skill yourself — it is deliberately *not* called from
`deploy.R`.) Re-run `document` if a fix touched roxygen.

## Step 4 — R CMD check

```bash
Rscript .claude/skills/deploy-bgnormr/deploy.R check
```

Read `$DEPLOY_OUT/check.log`. Fix every ERROR and WARNING, and NOTEs where
feasible; re-run until the summary reads `0 error(s) | 0 warning(s)`. This is
slow (2–3 min: it builds the package and runs the tests **and the vignette**).

## Step 5 — BiocCheck

```bash
Rscript .claude/skills/deploy-bgnormr/deploy.R bioccheck
```

Mirrors the CI call (`no-check-version-num = TRUE`). Fix all errors and
warnings; resolve notes where reasonable, else leave them (the current tree
carries 7 acceptable notes — line length, function length, `fnd` author role,
`dontrun` usage). See Gotchas for the two "errors" the helper marks
**environmental / non-blocking**.

## Step 6 — Version + NEWS

Bump the `z` in `Version:` (0.99.z convention — every push during Bioc review
bumps z) in [DESCRIPTION](DESCRIPTION), and add or extend the top entry in
[NEWS.md](NEWS.md) describing what this push changes.

## Step 7 — Commit & push

Stage, commit, and push to the branch confirmed in step 0. End the commit
message with:

```
Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
```

If the target is `main`, branch first per repo policy unless the user
explicitly asked to push straight to `main`.

## Gotchas

- **The `<pkg>.BiocCheck` output folder.** BiocCheck writes a `bgnormR.BiocCheck/`
  folder into the repo, then flags a leftover one on the *next* run
  (`checkBiocCheckOutputFolder` ERROR). `deploy.R` deletes it before and after
  each run so this never fires — but if you call `BiocCheck()` by hand,
  `unlink("bgnormR.BiocCheck", recursive = TRUE)` afterwards, and never commit it.
- **`checkSupportReg` "ERROR" is a network flake.** BiocCheck hits the Bioc
  support site to verify the maintainer's email; it fails with `HTTP 504` /
  "Unable to find your email" when offline or the site is slow. `deploy.R`
  classifies both this and the output-folder check as environmental and does
  **not** count them toward the gate. Don't chase them.
- **A broken vignette aborts `check` before findings exist.** The vignette is
  re-built inside `pkgbuild::build()`, *before* R CMD check runs, and that throws
  regardless of `error_on`. `deploy.R` catches it and reports
  `check: build aborted` — treat it as a blocking error and fix the vignette.
- **BiocCheck warnings block, not just errors.** For Bioconductor a WARNING is a
  release blocker; `deploy.R` counts check+BiocCheck warnings toward the gate too.
- **`document()` can change `NAMESPACE` unexpectedly.** Always diff it after step 1.

## Troubleshooting

- **`check: build aborted` → `could not find function "labs"` (or any bare
  ggplot2 verb) in a vignette chunk.** The vignette calls a ggplot2 function
  (`labs`, `theme`, `aes`, …) directly but only does `library(bgnormR)`. Those
  functions are *imported* into bgnormR's namespace, not re-exported, so a bare
  call can't resolve them. Fix: add `library(ggplot2)` to the vignette's setup
  chunk, or qualify the call (`ggplot2::labs(...)`). This is the exact failure
  the check phase caught while authoring this skill.
- **`there is no package called 'bgnormR'` from deploy.R.** Run it from the repo
  root; the helper uses `devtools`, which loads the source tree in place.
- **BiocCheck note count looks huge.** Each note lists every offending
  file:line; the helper's summary counts *distinct checks* via `res$getNum()`,
  which is the number that matters (currently 7).
