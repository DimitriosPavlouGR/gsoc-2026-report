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
    \max_{x} \quad & x_k \\
    \text{s.t.} \quad
    & Sx = 0, \\
    & l \leq x \leq u \quad\\
    \end{aligned}
$$

$$
    \begin{aligned}
    \min_{x} \quad & x_k \\
    \text{s.t.} \quad
    & Sx = 0, \\
    & l \leq x \leq u \quad\\
    \end{aligned}
$$

The same LPs also answer the second question at no additional cost. If $x_k^{\max}-x_k^{\min} \leq \varepsilon$, the reaction cannot vary within the feasible polytope and its therefore pinned.

To test whether the upper bound $u_k$ is essential, the algorithm moves it outwards by one and solves the maximization problem again:

$$
    \begin{aligned}
    \max_{x} \quad & x_k \\
    \text{s.t.} \quad
    & Sx = 0, \\
    & l_i \leq x_i \leq u_k \quad (i \neq k),\\
    & l_k \leq x_k \leq u_k+1
    \end{aligned}
$$

If the optimum does not increase, the original bound was not essential for the description of the polytope and can therefore be relaxed to $+\infty$. The lower bound can be tested simillary by relaxing $l_k$ to $l_k-1$ and solving a minimization LP. Thus, at most four LPs are required per reaction.

Fixing a reaction changes the feasible polytope, so the remaining bounds must be checked once again. A bound that appeared essential in an earlier iteration may become redundant after another constraint is removed. The algorithm therefore repeatedly goes over the reaction variables and terminates on the first complete pass where no bounds have been relaxed or a variable has been pinned.

```c++
MetabolicPolytope<Point> P;
...
ExhaustiveConfig config;
ExhaustiveSimplifier f(P, config);
MetabolicPolytope<Point> Ps = simplify(P, config).first;
```
##### Advantages and Limitations

The main advantage of `ExhaustiveSimplifier` is that it is robust to numerical issues. If a HiGHS LP fails for a particular bound, the algorithm can leave that bound unchanged and continue with the rest of the reactions. One problematic LP therefore does not necessarily stop the whole simplification process.

The main disadvantage is that this method can be slow. Each iteration can require up to four LP solves per pass, and several passes may be needed. Many LPs are solved using most of the bounds of the original model, with most of its variables and constraints. For a metabolic network with thousands of reactions and constraints, this leads to a very large number of large LP calls and makes the exhaustive approach expensive for larger models.

#### Clarkson's Simplification 
The second algorithm implemented tries to address the drawbacks of the exhaustive approach. Most bounds are not facets of the polytope, so solving an LP for every bound is wasteful. It keeps two sets of bounds, an essential set $I$, containing bounds that have been proved to define facets, and a set $J$, containing bounds whose status is still unknown. Initially $I$ is empty and $J$ contains all the constraints. All reactions start free, so the LPs only use the equality constraints $A_{eq} x = b_{eq} and the bounds I.

At each iteration, one bound is selected from $J$ and tested against the polytope described by $I$. The bound is relaxed and the corresponding reaction is optimized in the direction of the relaxed bound. If the optimum stays $\varepsilon$ within the original bound, then the bounds already in $I$ imply it, so the bound is redundant and can be removed. If the optimum lies outside the original polytope, the algorithm finds an essential constraint, i.e. a facet of
the polytope that blocks this direction. Starting from an interior point $z$, it follows the ray:

$$z+t(x^{\ast}-z)$$

to the optimum $x^{\ast}$ of the LP solved and finds the first bound that the ray hits. This bound must be a facet of the polytope, so it is added to $I$. This allows the algorithm to discover only the constraints that are actually needed to describe the polytope.

The interior point is found using the following LP:

$$\max\ y \quad \text{s.t.} \quad A_{eq}x = b_{eq}, \; b_l+y \leq x \leq b_u-y,\;0 \leq y \leq 1$$

Here, maximizing $y$ finds a point that is as far as possible from the finite bounds, providing an interior point for the ray shooting step.

Reaction pinning is also handled more efficiently. Instead of repeatedly checking every reaction, the implementation maintains a worklist and only requeues reactions that share an equality with the row corresponding to a reaction that was just fixed.

```c++
MetabolicPolytope<Point> P;
...
ClarksonConfig config;
ClarksonSimplifier f(P, config);
MetabolicPolytope<Point> Ps = simplify(P, config).first;
```

##### Advantages and Limitations

The main advantage is that LPs are solved against the much smaller set $I$ rather than against every bound in the original model. Since most bounds are usually redundant, this can greatly reduce the number of constraints and LP calls compared with the exhaustive approach.

The main disadvantage is that the method is more complex and relies on additional geometric and numerical operations, such as finding an interior point and determining which facet a ray intersects. This makes the implementation much more sensitive to numerical issues than the exhaustive approach.

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

### Replacing LPSolve

The lpsolve calls weren't confined to `include/lp_oracles`, they were reached indirectly by the polytope classes, the random walks, and then volume algorithms, so every oracle had to be traced carefully through its callers before I could touch it. In addition HiGHS.cmake had to follow the pattern VolEsti's other dependencies already used, and before converting anything I had to settle on the new signatures and implementations of the oracles.

### Implementing Exhaustive & Clarkson 

A bug in either algorithm doesn't crash the program, but it can quietly hand back incorrect polytopes, so I couldn't take that the code worked properly for granted even if it seemed correct. That meant implementing additional oracles
for testing containment and equallity between different metabolic polytopes. The difficult part
was telling apart three different types of failures, a real bug, a numerical decision made too close to the solver's own tolerance, and an LP that failed because the model was badly scaled, which is what motivated the scaling stage. Comparing the two algorithms fairly also meant judging them by the region they describe rather by raw counts of bounds relaxed or reactions pinned, since the same region can be reached through different operations with different counts.

### Contributing to open source

The code had to read as part of VolEsti, not bolted onto it, so I tried to match its conventions where they existed, e.g. (no namespaces) rather than doing things my way, and wrote plenty of tests to make sure everything worked properly.

## Acknowledgments

I would like to thank my mentors Vissarion Fisikopoulos and Apostolos Chalkis for giving the opportunity to work on VolEsti and for their constant guidance and support throughout GSoC 2026.
