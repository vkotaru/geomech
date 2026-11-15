# Sums over n bodies

A chain of $n$ links, a cable modeled by $n$ lumped masses, or a team of $n$ quadrotors
has one equation of motion per body, and the equations of the bodies share one form.
`geomech.indexed` derives that form once: the equation of the $i$-th body, with $i$ and
$n$ symbolic, written with sums whose bounds depend on $i$.

## Writing a model

A `Family` is a set of symbols $q_1, \dots, q_n$ of one kind; `Sum(j, lower, upper, body)`
is $\sum_{j=\text{lower}}^{\text{upper}}$ body. A bound is `1` (or another integer), the
number of bodies `n`, an index, or `i - 1` / `i + 1` for a strict bound. The n-link
spherical chain, with $x_i = \sum_{j \le i} l_j q_j$:

```python
from geomech import *
from geomech.indexed import Family, IndexedVariables, Sum, compute_indexed_eom, indices, n

i, j = indices("i j")
g = Scalar("g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
q = Family(S2, "q")                              # q_1, ..., q_n on S^2
m = Family(Scalar, "m", attr=["Constant"])
l = Family(Scalar, "l", attr=["Constant"])
f = Family(Vector, "f")

def x(k):                                        # position of mass k
    return Sum(j, 1, k, q[j] * l[j])

v = x(i).t_diff()
L = Sum(i, 1, n, 0.5 * m[i] * v.dot(v)) - Sum(i, 1, n, m[i] * g * x(i).dot(e3))
dW = Sum(i, 1, n, Dot(q[i].get_variation_vector(), f[i]))

eom = compute_indexed_eom(L, dW, IndexedVariables(families=[q]))
print(eom["q"])
```

A sum supports `+`, `-`, `*` (with expressions, numbers and other sums), `.dot`,
`.cross`, `.t_diff()` and `.T()`; the product of two sums is a double sum over two
indices. `IndexedVariables` lists the families of configuration variables and, as
`scalars`, `vectors`, `matrices`, the ordinary ones (the quadrotor of a flexible cable).

## The result

`eom["q"]` is the equation of $q_i$:

$$
f_i \;-\; \sum_{j=i}^{n} m_j\, g\, l_i\, (q_i \times e_3)
\;+\; \sum_{j=1}^{n}\, \sum_{k=\max(i,j)}^{n} m_k\, l_i\, l_j \big(\cdots\big) = 0 .
$$

The inner sum $\sum_{k \ge \max(i,j)} m_k\, l_i l_j$ is the inertia $M_{ij}$ of Lee et
al. §5.3: the mass of every link below both $i$ and $j$. The equation of an ordinary
variable keeps its sums over all bodies; for the flexible cable, the quadrotor's
translation is driven by $m_Q$ and by $\sum_j m_j$.

`eq.to_latex()` prints the equation with nested `\sum`; `eq.instantiate(size, at=k)`
expands it for $n$ = `size` at body $k$, and `eom.instantiate(size)` gives the
dictionary `compute_eom` would return for that $n$, which `check_eom`, `to_standard_form`,
`linearize` and the code generator take as usual.

## How it works

A sum is not a node of the expression tree. It is kept as a list of *terms*, each a
summand (an ordinary expression whose leaves carry the index in their name, `q_{j}`)
with its bound indices and the order constraints between them ($1 \le j \le i \le n$).
`compute_indexed_eom` takes each term through the stages of `compute_eom` on its own
(variation, manifold rules, simplification, integration by parts, extraction), since
the variation and the time derivative commute with a sum of fixed length. The
coefficient of $\xi_j$ under a bound index $j$ becomes a term of the equation of body
$i$ by the substitution $j := i$, and the constraints on $j$ become constraints on the
other indices: $\sum_{k} [\,i \le k\,] \cdots = \sum_{k \ge i} \cdots$. Terms that differ
only in the names of their bound indices are merged, and each summand is simplified.

## Verification

For every model in `geomech.indexed.systems` (the spherical and rigid chains, the
quadrotor with a flexible cable, and $n$ quadrotors carrying a point mass or a rigid
body), the tests expand the generic equation for $n = 1, \dots, 6$ and run
`check_eom` against the Lagrangian expanded for the same $n$; removing any one derived
term makes the check fail. For the spherical chain, the generic equation is also
compared with Lee et al. eq. 5.17 at every $n$ tested.

## Limits

All families share one number of members $n$, and a member is indexed by an index or an
integer (`q[i]`, `q[3]`), not by `i - 1`. A model whose virtual work couples neighbors
through `i - 1` (joint torques shared by two links) is written with body moments instead.
`to_standard_form`, `linearize` and the code generator work on an instantiated model
(a given $n$).
