# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

ThreeBodyProblem.jl is a registered Julia package (General registry) for astrodynamics in the Circular Restricted
Three-Body Problem (CR3BP), with supporting two-body, Clohessy-Wiltshire, bicircular (BCP) and elliptic (ER3BP)
models. Docs: https://jared711.github.io/ThreeBodyProblem.jl. Worked examples live in a separate repo,
[ThreeBodyProblemExamples.jl](https://github.com/jared711/ThreeBodyProblemExamples.jl).

## Commands

Julia is managed with juliaup (1.12 locally). This checkout is `Pkg.develop`ed into the default Julia env, so
`using ThreeBodyProblem` anywhere loads this working tree.

```bash
julia --project -e 'using Pkg; Pkg.test()'      # full suite (test/runtests.jl, one @safetestset per src file)
julia --project test/util_tests.jl              # one test file; files are self-contained (using ThreeBodyProblem, Test)
julia --project=docs -e 'using Pkg; Pkg.instantiate()'   # first time only (docs/Project.toml uses [sources] path="..")
julia --project=docs docs/make.jl               # build docs into docs/build; also runs the jldoctests
```

The docs build is strict (Documenter 1.x): every docstring must appear in an `@docs`/`@autodocs` block on some
page in `docs/src/`, every name listed there must exist, and every ```` ```jldoctest ```` output must match
exactly. When adding or renaming a documented function, also update the matching `docs/src/*.md` page.
Plots is a test-only dependency (`[extras]`/`[targets]`); the package itself depends only on RecipesBase.

`test/special_tests.jl` is commented out in `runtests.jl`: it calls functions (`rich3`, `differential_corrector`,
`invariant_manifolds`) that are not in the package.

## Architecture

`src/ThreeBodyProblem.jl` `include`s each file and owns the explicit `export` list. A new public function must be
added there too. The package has no integrator dependency (a work-in-progress Julia port of a MATLAB RKF7(8)
integrator, `ode78e`, sits in the gitignored `example/integrators.jl`): the dynamics functions are written in
DifferentialEquations.jl's `f(u, p, t)` / `f!(du, u, p, t)` signature, and users integrate them with OrdinaryDiffEq.

- `parameters.jl`: `Body` (m, R, a, T, name, color) and the system structs `System` (CR3BP), `BicircularSystem`,
  `EllipticalSystem`. A `System` stores both ASCII and Unicode copies of the same fields (`mu`/`μ`, `mu1`/`μ₁`,
  `R1`/`R₁`) plus normalizing units `RUNIT, TUNIT, VUNIT, AUNIT, MUNIT`. Body constants (`EARTH`, `MOON`,
  `SATURN`, `ENCELADUS`...) and constructors (`earth_moon()`, `saturn_enceladus()`...) are defined here. Units are
  km, s, kg.
- `dynamics.jl`, `util.jl`, `frames.jl`: most functions come in several dispatch variants with the same name:
  - the third argument is a `μ::Real` (nondimensional, rotating frame), a `sys::System`, or a `p::Array = [μ₁, μ₂, d]`
    (dimensional, in km and km³/s²);
  - an in-place `foo!` version writes into its first argument and returns `nothing` (tests check `isnothing(...)`).
  Keep all variants consistent when changing one. State vectors are `rv = [x,y,z,ẋ,ẏ,ż]`. STM propagation uses
  `w = [rv; vec(Φ)]` (length 42), with `CR3BPjac` as the variational matrix.
- `plot.jl`: RecipesBase `@recipe`s for `System`, `Body` and trajectories (`Vector{Vector{<:Real}}`), so
  `plot(sys)` works once the user loads Plots. Only the helper geometry (`circle`, `sphere`, `torus`, `seczoom`) is
  exported.
- `util.jl` also holds Lagrange points (`computeL1`..`computeLpts`), Jacobi constant (`computeC`, `rC2v`),
  `eig` (eigen-decomposition ordered for monodromy matrices), `stability_index`, rotation matrices and
  time/angle helpers.

## Downstream user: nrho-lagrangian-tori

`~/Projects/nrho-lagrangian-tori` (QPO / invariant-torus families around Saturn–Enceladus NRHOs) depends on the
**registered** v1.0.3, not on this checkout. Before releasing, check that it still works with these APIs:
`System(ENCELADUS, SATURN)` and the **full 20-argument positional `System` constructor**
(`src/julia/main.jl:140`). Reordering or adding `System` fields breaks it. It also uses `plot!(sys; sec, Lpts,
lims, center)`, `eig`, `stability_index`, `ThreeBodyProblem.G` and `ThreeBodyProblem.rotz`.

## CI and releases

GitHub Actions: `CI.yml` (tests on Julia LTS and current release, plus Windows/macOS; coverage to Codecov),
`Documentation.yml` (builds docs and deploys to `gh-pages` via `deploydocs`), and `TagBot.yml`. Releases are made
by bumping `version` in `Project.toml` and registering it with the General registry (Registrator), after which
TagBot tags the release and the tag build publishes the `stable` docs.

`example/`, `src/fig/`, `Manifest.toml` and `docs/Manifest.toml` are gitignored local files, not part of the package.
