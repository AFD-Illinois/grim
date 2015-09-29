---
layout: page
title: Motivation
---

#### Why $$\mathtt{grim}$$?
---
Conservative, shock-capturing codes solve for a chosen set of
_primitive_ variables $$P$$ , such as density, pressure, and so on, by evolving
the governing equations in their _conservative_ form
\begin{align}
\frac{\partial U}{\partial t} + \nabla\cdot F & = S
\end{align}
where $$U\equiv U(P)$$ are  _conserved_ variables, $$F\equiv F(P)$$ are spatial
_fluxes_ and $$S \equiv S(P)$$ are _sources_. The time evolution is performed
over the conserved variables $$U$$, i.e. $$U^n\rightarrow U^{n+1}$$, where $$n,
n+1$$ indicate the discrete time levels. Solving for the primitive variables
$$P^{n+1}$$ from the time evolved conserved variables $$U^{n+1}(P^{n+1})$$ at
the end of a discrete time step is a nonlinear root finding problem. The
primitive variable recovery algorithms used in existing codes only work for
perfect fluids. They cannot be used for the governing equations that include
dissipative processes that arise when one deviates from local thermodynamic
equilibirum[^EMHD_model_paper].  Hence the need arose to write a physics
agnostic code that can be used to explore a variety of models, independent of
the sophistication of the model.

[^EMHD_model_paper]: Chandra, Gammie, Foucart, Quataert (2015), [arXiv:1508.00878](http://arxiv.org/abs/1508.00878)

#### How does $$\mathtt{grim}$$ do it?
---
In $$\mathtt{grim}$$, the governing equations of a particular model are written
down in their conserved form and recast as a massive coupled nonlinear root
finding problem. Writing down the equations of the model explicitly in terms of
primitive variables to be solved for, $$P^{n+1}$$,
\begin{align}
R(P^{n+1}) = \frac{\partial U}{\partial t} + \nabla\cdot F - S
\end{align}
where $$R(P^{n+1})$$ are the _residuals_. We wish to find those $$P^{n+1}$$ for
which $$R(P^{n+1}) = 0$$. We do so using the _Newton-Krylov_ algorithm, which
iteratively solves for $$R(P^{n+1})\lt tol$$, where $$tol$$ is a chosen error
tolerance. The algorithm numerically assembles the Jacobian during the iteration
process using only the residuals $$R(P^{n+1})$$ as inputs. This automatic
Jacobian assembly enables the inclusion of any physical model, including those
in which the source terms have spatio-temporal derivatives. 
$$\mathtt{grim}$$ is built on top of the
$$\mathtt{PETSc}$$[^PETSc_webpage] framework and uses the Newton-Krylov implementation in the
$$\mathtt{snes}$$ module. The spatio-temporal discretizations are performed over
a structured grid constructed using the $$\mathtt{DMDA}$$ module. Using the
$$\mathtt{DMDA}$$ module to construct the grid gives $$\mathtt{snes}$$
information about the width of the reconstruction stencils being used, and the
connectivity between neighbouring spatial elements. With this information,
$$\mathtt{PETSc}$$ can use graph coloring and efficiently assemble the Jacobian.
This allows us to use both _explicit_ $$\nabla\cdot F \equiv (\nabla\cdot F)
(P^{n})$$, as well as _implicit_ $$\nabla \cdot F \equiv (\nabla \cdot
F)(P^{n+1}, P^{n})$$ time stepping, and merge the code needed to do both into a
single residual assembly routine. For details, see the
[algorithms](../algorithms) page.

[^PETSc_webpage]: [Portable, Extensible Toolkit for Scientific Computation](http://www.mcs.anl.gov/petsc/)
