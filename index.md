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

**Simplification for metabolic networks:** A metabolic network is described by a stoichiometric matrix $S$ and per reaction flux bounds, giving the polytope: 

$$P=\{x \in \mathbb{R}^d : Sx = 0,\; b_l \leq x \leq b_u\}$$

whose points are the steady-state flux distributions of the network. Genome scale models make this polytope awkward to hand to sampling/volume approximation algorithms directly, since most of the flux bounds turn out to be redundant or some reactions turn out to be pinned to a single value once the stoichiometric matrix is taken into account. I added a `MetabolicPolytope` representation with a parser for BiGG models, along with a pipeline that prepares it for sampling/volume approximation. A scaling stage is used internally to bring the bounds and coefficients to comparable magnitudes. Simplification then relaxes the redundant bounds and moves the pinned reactions into $S$, either by the exhaustive or Clarkson's algorithm. What remains is transformed into the full dimensional H-polytope that VolEsti's existing volume and sampling algorithms accept.

## Pull Requests

| PR | Description | Status |
|----|-------------|--------|
| [#493](https://github.com/GeomScale/volesti/pull/493) | Metabolic polytope representation, and exhaustive simplification | Closed|
| [#499](https://github.com/GeomScale/volesti/pull/499) | Finished simplification, added Clarkson, BiGG parser, and exhaustive simplification | Merged |
| [#XXX](https://github.com/GeomScale/volesti/pull/XXX) | Added scaling, and HiGHS replacing lpsolve across the LP oracles | Open |

## Summary of Contributions
### Replacing lpsolve with HiGHS
Every sampling and volume algorithm in VolEsti sits on top of a handful of linear programming oracles located in `include/lp_oracles`, which locate an interior point, test membership, and compute where a ray leaves a body, among other things. All of them called lpsolve, a solver that is no longer maintained and that loses both speed and accuracy in high dimensions.

**Updating the build system:** lpsolve was removed and `HiGHS.cmake` was added following the pattern the other dependencies already used, so HiGHS is fetched and built on demand and no manual installation is needed, with the examples and tests linking it through `highs::highs` target.

**Updating the oracles:** Each oracle builds a model, solves it, and reads back either the objective or the solution, so the migration was a rewrite of the model construction in HiGHS's API while preserving the behaviour of each one.

| Oracle | Computes |
|--------|----------|
|`compute_chebychev_ball` | the largest ball inscribed in an H-polytope |
|`point_in_intersection` |a point common to two V-polytopes |
|`memLP_Vpoly`, `memLP_Zonotope` | whether a point belongs to a V-polytope or a zonotope |
|`intersect_line_Vpoly` | the point where a ray hits a V-polytope |
|`intersect_double_line_Vpoly`, `intersect_line_zono` | both intersections of a line with a V-polytope or a zonotope |
|`is_contained_in`, `are_equal` | whether one metabolic polytope contains, or equals another|

**Configuring the solver:** Every oracle now takes an optional `LPOracleOptions`, a callable applied to the `Highs` instance before the model is built and solved, so a caller can set a time limit, a tolerance, choose a solver, or pick any other configuration without the oracle having to expose each option itself.

```cpp
LPOracleOptions opts = [](Highs& highs) {
    highs.setOptionValue("time_limit", 5.0);
    highs.setOptionValue("solver", "simplex");
};

auto ball = compute_chebychev_ball<NT, Point>(A, b, opts);
```

**Error handling:** The oracles returned their result directly, without additional information on whether the solve had succeeded. Each one now returns an `LPOracleResult<T>` struct, which carries the computed value together with a `solved` flag, so the value is only read once the solve is confirmed to be successful.
```cpp
auto res = point_in_intersection<VT>(V1, V2, direction);

if (!res.solved) {
    // the LP failed
} else if (res.value.second) {
    // the polytopes do not intersect
} else {
    // res.value.first is a point in the intersection
}
```
### Metabolic Polytope Preprocessing
#### Metabolic Polytope 
A new class called `MetabolicPolytope` was added that stores a metabolic network in the form it naturally comes in, box bounds together with equality constraints, rather than as a single system of inequalities:

$$P=\{x \in \mathbb{R}^d : Sx = 0,\; b_l \leq x \leq b_u\}$$

The bounds are held as two vectors and the equalities as a matrix that is sparse and row-major by default, which suits the stoichiometric matrix of a genome scale model. The class is templated on the point type and on the matrix type, so the storage can be changed where a different one fits better.

```c++
template
<
    typename Point,
    typename MT_Type = Eigen::SparseMatrix<typename Point::FT, Eigen::RowMajor>
>
class MetabolicPolytope;
```
A metabolic network enters this form with $A_{eq} = S$ and $b_{eq} = 0$, the steady state condition. The right hand side is kept as a vector rather than fixed to zero so that the simplification algorithm can add discovered equalities to $A_{eq}$ when simplifying the metabolic network.

**A parser for BiGG models:** Models from [BiGG](http://bigg.ucsd.edu/) are distributed as SMBL, MAT and JSON files. So a parser was added that reads a BiGG model in JSON format and turns it into a `MetabolicPolytope` using [nlohmann/json](https://github.com/nlohmann/json).

```c++
auto P = parse_from_json<Point>("e_coli_core.json");
```
#### Simplification Overview
Once a model is loaded it cannot be handed to VolEsti's sampling and volume approximation implementations directly. The equalities present in the stoichiometric matrix $A_{eq}$ coffine it to an affine subspace of $\mathbb{R}^d$, so it has no interior and zero volume in the ambient space, while VolEsti's implementations expect a full dimensional H-polytope. The models also often happen to be far more complex than they need to be, the steady state condition makes most of the inequalities redundant, which means that they are not facets of $P$ and only add rows to the system that the algorithms will keep checking. Furthermore, some reactions are constrained so tightly by the network that they cannot vary at all, and so are not free variables, but hidden equalities, reducing the actual dimension of the polytope even further.

Both problems are addressed before transforming the polytope to a full dimensional polytope. Simplification relaxes the redundant bounds and moves the pinned reactions into $A_{eq}$, which shrinks the description of the polytope. Two algorithms have been implemented to simplify the polytopes, `ExhaustiveSimplifier`, and `ClarksonSimplifier`.

#### Exhaustive Simplification 
The first algorithm answers both questions by solving up to four LPs for every single reaction. For a reaction $(k)$, it first determines its maximum and minimum feasible flux by solving:

$$
    \begin{aligned}
    \max_{x},\quad \min_{x} \quad & x_k \\
    \text{s.t.} \quad
    & Sx = 0, \\
    & l_i \leq x_i \leq u_l \quad\\
    \end{aligned}
$$

The same LPs also answer the second question at no additional cost. If $x_k^{\max}-x_k^{\min} \leq \varepsilon$, the reaction cannot vary within the feasible polytope and its therefore pinned.

To test whether the upper bound $u_k$ is essential, the algorithm moves it outwards by one and solves the maximization problem again:

$$
    \begin{aligned}
    \max_{x} \quad & x_k \\
    \text{s.t.} \quad
    & Sx = 0, \\
    & l_i \leq x_i \leq u_l \quad (i \neq k),\\
    & l_k \leq x_k \leq u_k+1
    \end{aligned}
$$

The maximization and minimization LPs also answer the second question at no extra cost. If the two optimums are close within $\epsilon$, the reaction cannot vary and is pinned.

Relaxing a bound or fixing a reaction changes the polytope, therefore the remaining bounds are checked again, and can expose bounds that looked essential in previous iteration. The algorithm therefore sweeps the reactions repeatedly and stops on the first pass that changes nothing.

```c++
using namespace exhaustive_simplification;
Config config;
config.fix_dimensions = true;
config.verbosity = VerbosityLevel::Summary;

auto res = simplify(P, config);
```
The drawback of this implementation is that every LP is solved against the whole model, and up to four per reaction LPs per pass are needed, which on a genome scale network is a lot of LP calls over a system with thousands of columns and rows.

#### Clarkson's Simplification 
The second algorithm implemented tries to address the drawbacks of the first implementation, which is that every LP is solved against every bound even though almost none of them describe the facet. It keeps two sets: an essential set $I$ of bounds proved to be essential, i.e. facets, initially empty, and a set $J$ of bounds whose status is still unknown, initially holding every finite bound of the input. Every reaction of the model (column) starts free, so the LPs are solved against $A_{eq}$ and $I$ alone, and $I$ only ever contains genuine facets.

Each iteration draws a bound from $J$ and asks the same question the exhaustive method asks, but against $I$ instead of the full polytope, and with a single LP call. The bound is applied to the model relaxed outwards and the reaction is optimized in its own direction. If the optimum lands within the original bound, then the facets already in $I$ already imply it, and the bound is deemed redundant and dropped. If the optimum $x^\ast$ lands outside the polytope, then the segment from an interior point $z$ of the polytope towards it must leave through a facet of the polytope. The ray $z+t(x^*-z)$ is tested against every box bound and the first facet hit is added to $I$.

The interior point is found using the following LP:

$$\max\ y \quad \text{s.t.} \quad A_{eq}x = b_{eq}, \; b_l+y \leq x \leq b_u-y,\;0 \leq y \leq 1$$

Reaction pinning is also done differently here as well. Instead of going over all the reactions repeatedly, the implementation keeps a worklist and requeues only the reactions sharing an equality with the row that was just fixed.

```c++
using namespace clarkson_simplification;
Config config;
config.fix_dimensions = true;
config.verbosity = VerbosityLevel::Summary;

auto res = simplify(P, config);
```

#### Scaling 
The numbers in a genome scale model do not share a magnitude. Flux bounds range from very large to very tiny values, and the stoichiometric coefficients have a spread of their own. A scaling stage was added to bring everything onto a comparable ground before the LPs are solved, which makes the tolerances meaningful.

A scaling is two strictly positive vectors, one factor per reaction and one per metabolite:

$$A^\ast_{ij} = \frac{A_{ij}}{r_i c_j}, \quad b^\ast_{eq,i} = \frac{b_{eq,i}}{r_i},\quad x^\ast _j = c_j x_j,\quad b^\ast_{l,j}=c_j b_{l,j},\quad b^\ast_{u,j}=c_j b_{u,j}$$

This is a change of variables, so the structure of the polytope is mostly untouched and only the numerical values move. The choice of factors is left to a policy, a callable that takes the polytope and fills in the scaling vectors $row,col$.

```c++
Scaling<Point> s;
auto Ps = scale(P, s, MaxBoundScaling{});

auto x = scale_point<Point>(x_scaled, s, false); // back to original coordinates
```
Two scalings have also been implemented, `MaxBoundScaling`, which divides each reaction by its largest finite bound and each metabolite row by its largest coefficient, sending both to $\pm 1$. `GMScaling` follows the geometric mean approach of `gmscale`, alternating column and row passes until the passes stop improving the factors.

## Challenges
**Working in a large codebase:** VolEsti is a template heavy library, and the LP oracles sit underneath almost everything in it. Understanding how a change would affect the other components of the library meant reading well beyond the files I was editing, since the polytope classes, the random walks and the volume algorithms all reach the oracles indirectly.

**Replacing lpsolve:** The lpsolve calls were spread across the oracles with no abstraction, so the migration had to preserve the behaviour of each one exactly while changing the solver underneath.

## Acknowledgments

I would like to thank my mentors for their guidance and support throughout GSoC 2026.
