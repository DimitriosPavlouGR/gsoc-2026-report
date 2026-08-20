---
title: Overview
---

# Metabolic network simplification and a new LP backend for VolEsti

**Contributor:** Dimitrios Pavlou<br>
**Mentors:** Vissarion Fisikopoulos, Apostolos Chalkis<br>
**Organization:** [Geomsacle](https://geomscale.github.io/about/)

[VolEsti](https://github.com/GeomScale/volesti) is a C++ library for volume approximation and sampling of convex bodies (e.g. polytopes) with an R interface.

## Project Overview
During Google Summer of Code 2026, I made two related contributions to VolEsti.

**Replacing lpsolve with HiGHS:** lpsolve is obsolete, slow, and loses accuracy in high dimensions, which matters because every oracle built on it feeds a volume or sampling algorithm. Every LP oracle in the library, e.g. Chebychev ball computation, V-polytope and zonotope membership, ray and line intersection, and V-polytope intersection now runs on HiGHS, behind an interface that lets callers configure the solver and that reports solver failure and infeasibility separately.

**Simplification for metabolic networks:** A metabolic network is described by a stoichiometric matrix `S` and per reaction flux bounds, giving the polytope: 

$$P=\{x \in \mathbb{R}^d : Sx = 0,\; b_l \leq x \leq b_u\}$$

whose points are the steady-state flux distributions of the network. Genome scale models make this polytope awkward to hand to sampling/volume approximation algorithms directly, since most of the flux bounds turn out to be redundant once the stoichiometric matrix is taken into account, some reactions are pinned to a single value. I added a `MetabolicPolytope` representation with a parser for BiGG models, along with a pipeline that prepares it for sampling/volume approximation. A scaling stage is used internally to bring the bounds and coefficients to comparable magnitudes. Simplification then relaxes the redundant bounds and moves the pinned reactions into `S`, either by the exhaustive or Clarkson's algorithm. What remains is transformed into the full dimensional H-polytope that VolEsti's existing volume and sampling algorithms accept.

## Pull Requests

| PR | Description | Status |
|----|-------------|--------|
| [#493](https://github.com/GeomScale/volesti/pull/493) | Metabolic polytope representation, and exhaustive simplification | Closed|
| [#499](https://github.com/GeomScale/volesti/pull/499) | Finished simplification, added Clarkson, BiGG parser, and exhaustive simplification | Merged |
| [#XXX](https://github.com/GeomScale/volesti/pull/XXX) | Added scaling, and HiGHS replacing lpsolve across the LP oracles | Open |

## Summary of Contributions
### Replacing lpsolve with HiGHS
Every sampling and volume algorithm in VolEsti sits on top of a handful of linear programming oracles located in `include/lp_oracles`, which locate an interior point, test membership, and compute where a ray leaves a body, among others. All of them called lpsolve, a solver that is no longer maintained and that loses both speed and accuracy in high dimensions.

**Updating the build system:** lpsolve was removed and `HiGHS.cmake` was added following the pattern the other dependencies already used, so HiGHS is fetched and built on demand and no manual installation is needed, with the examples and tests linking it through `highs::highs` target.

**Updating the oracles:** Each oracle builds a model, solves it, and reads back either the objective or the solution, so the migration was a rewrite of the model construction in HiGHS's API while preserving the behaviour of each one.

| Oracle | Computes |
|`ComputeChebychevBall` | the largest ball inscribed in an H-polytope |
|`memLP_Vpoly`, `memLP_Zonotope` | whether a polytope belongs to a V-polytope or a zonotope |
|`intersect_line_Vpoly` | the point where a ray hits a V-polytope |
|`intersect_double_line_Vpoly`, `intersect_line_zono` | both intersection of a line with a V-polytope or a zonotope |
|`PointInIntersection` |a point common to two V-polytopes |
<br><br>
**Configuring the solver:** Every oracle now takes an optional `LPOracleOptions`, a callable applied to the `Highs` instance before the model is built and solved, so a caller can set a time limit, a tolerance, choose a solver, or pick any other configuration without the oracle having to expose each option itself.

```cpp
LPOracleOptions opts = [](Highs& highs) {
    highs.setOptionValue("time_limit", 5.0);
    highs.setOptionValue("solver", "simplex");
};
```

**Error handling:** The oracles returned their result directly, without additional information on whether the solve had succeeded. Each one now returns a named struct carrying the result together with that information, and where the distinction is meaningful an infeasible LP is reported as a fact about the geometry rather than as a failiure.

```cpp
auto res = PointInIntersection<VT>(V1, V2, direction);

if (!res.success) {
    // the LP failed
} else if (res.is_empty) {
    // the polytopes do not intersect
}
```
### Simplification


## Challenges
**Working in a large codebase:** VolEsti is a template heavy library, and the LP oracles sit underneath almost everything in it. Understanding how a change would affect the other components of the library meant reading well beyond the files I was editing, since the polytope classes, the random walks and the volume algorithms all reach the oracles indirectly.

**Replacing lpsolve:** The lpsolve calls were spread across the oracles with no abstraction, so the migration had to preserve the behaviour of each one exactly while changing the solver underneath.

## Acknowledgments

I would like to thank my mentors for their guidance and support throughout GSoC 2026.
