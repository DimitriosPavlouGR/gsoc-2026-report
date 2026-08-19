---
title: Overview
---

# Metabolic network simplification and a new LP backend for VolEsti

[VolEsti](https://github.com/GeomScale/volesti) is a C++ library for volume approximation and sampling of convex bodies (e.g. polytopes) with an R interface.

## Project Overview
During Google Summer of Code 2026, I made two related contributions to VolEsti.

**Replacing lpsolve with HiGHS:** lpsolve is obsolete, slow, and loses accuracy in high dimensions, which matters because every oracle built on it feeds a volume or sampling algorithm. Every LP oracle in the library, e.g. Chebychev ball computation, V-polytope and zonotope membership, ray and line intersection, and V-polytope intersection now runs on HiGHS, behind an interface that lets callers configure the solver and that reports solver failure and infeasibility separately.

**Simplification for metabolic networks:** A metabolic network is described by a stoichiometric matrix `S` and per reaction flux bounds, giving the polytope: 

$$P=\{x \in \mathbb{R}^d : Sx = 0,\; b_l \leq x \leq b_u\}$$

whose points are the steady-state flux distributions of the network. Genome scale models make this polytope awkward to hand to sampling/volume approximation algorithms directly, since most of the flux bounds turn out to be redundant once the stoichiometric matrix is taken into account, some reactions are pinned to a single value. I added a `MetabolicPolytope` representation with a parser for BiGG models, along with a pipeline that prepares it for sampling/volume approximation. A scaling stage is used internally to bring the bounds and coefficients to comparable magnitudes. Simplification then relaxes the redundant bounds and moves the pinned reactions into `S`, either by the exhaustive or Clarkson's algorithm. What remains is transformed into the full dimensional H-polytope that VolEsti's existing volume and sampling algorithms accept.

