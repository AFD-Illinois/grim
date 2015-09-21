---
layout: page
title: Motivation
---

#### Why $$\mathtt{grim}$$?
---
Conservative, shock-capturing codes solve equations that describe a chosen set
of _primitive_ variables $$P$$ , such as density, pressure, and so on, in the
following
_conservative_ form
\begin{align}
\frac{\partial U}{\partial t} + \nabla\cdot F & = S
\end{align}
where $$U$$ are  _conserved_ variables, $$F$$ are spatial _fluxes_ and $$S$$ are
_sources_. The time evolution is described in terms of the conserved variables
$$U$$. Recovering the primitive variables from the time evolved conserved
variables at the end of a discrete time step is a nonlinear root finding
problem. The recovery algorithms used in existing codes only work for perfect
fluids. They cannot be used for the governing equations that include dissipative
processes that arise when one deviates from local thermodynamic equilibirum.
Hence the need arose to write a physics agnostic code that can be used to
explore a variety of models, independent of the sophistication of the model.

#### How does $$\mathtt{grim}$$ do it?
---
In $$\mathtt{grim}$$, the governing equations of a particular model are written
down in their conserved form and recast as a massive coupled nonlinear root
finding problem. Writing $$U \equiv U(P)$$, $$F \equiv F(P)$$ and $$S \equiv
S(P)$$,
\begin{align}
R(P) = \frac{\partial U(P)}{\partial t} + \nabla\cdot F(P) - S(P)
\end{align}
where $$R(P)$$ are the _residuals_. We wish to solve for $$R(P) = 0$$. We do so
using the _Newton-Krylov_ algorithm, which iteratively solves for $$R(P)\lt
tol$$, where $$tol$$ is a chosen error tolerance. The algorithm automatically
assembles the Jacobian using colored finite differences during the iteration
process using only the residuals $$R(P)$$ as inputs. This automatic Jacobian
assembly enables the inclusion of any physical model, as well as _explicit_ and
_implicit_ time stepping schemes.
