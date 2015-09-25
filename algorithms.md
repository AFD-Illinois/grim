---
layout: page
title: Algorithms
---

$$\mathtt{grim}$$ uses a second order finite volume method to solve partial
differential equations in the following conservative form
\begin{align}
\frac{\partial U}{\partial t} + \nabla\cdot F & = S
\end{align}
in a spatial volume $$d$$ with the boundary $$\partial d$$. The formulation
starts by discretizaing the domain $$d$$, and is described below.

#### Grid Generation
---
The spatial volume $$d$$, is discretized by
transforming into a set of chosen computational coordinates $$X^\mu$$, in which
the coordinate axes are aligned to the transformed boundary $$\partial D$$. The
domain in these coordinates is then discretized using a structured uniform
hexahedral mesh.
![grid](../grid.png){: style="max-width: 500px; height: auto;"}

A common case involves discretizing spherical domains. This is used for example,
to study accretion flows. To construct a spherical mesh, a computational
coordinate system $$X^\mu$$ is chosen, which relates to the $$x^1$$ coordinates
by
$$\begin{align}
X^1 = \log(r) \\
\theta = \pi X^2 + \frac{1 - h}{2} \sin(2 \pi X^2)
\end{align}$$<br>
where $$r = \sqrt{x_1^2 + x_2^2}$$ and $$\theta = \tan^{-1}({x^2/x^1})$$. With
this coordinate mapping, a uniform grid in the $$X^\mu$$ coordinates leads to an
exponential packing of the grid zones at the inner boundary and a concentration
of the grid zones at the midplane which is controlled by the parameter $$h$$.
The discretized domain is illustrated below for $$h = 0.3$$.

![spherical_grid](../spherical_grid_generation.png){: style="max-width: 650px; height: auto;"}

#### Finite Volume Method
---
Muliplying the conservative equations by the volume element of a discrete zone
in the $$X^\mu$$ coordinate system, $$\Delta v = dX^1dX^2$$, and using the
divergence theorem, 
\begin{align}
\partial_t \bar{U} + \frac{\bar{F}^1_{right} - \bar{F}^1_{left}}{\Delta X^1} + \frac{\bar{F}^2_{top} - \bar{F}^2_{bottom}}{\Delta X^2} = \bar{S}
\end{align}
where $$\bar{U} = \int U \Delta v/\int \Delta v$$, $$\bar{S} = \int S \Delta v/\int
\Delta v$$ and $$\bar{F}^1 = \int F^1 dX^2/\int dX^2$$, $$\bar{F}^2 = \int F^2
dX^1/\int dX^1$$. The locations $$right$$, $$left$$, $$top$$ and $$bottom$$ are
illustrated below for an individual grid zone.

![gridzone](../gridzone.png){: style="max-width: 500px; height: auto;"}


#### Reconstruction
---
![reconstruction](../reconstruction.png){: style="max-width: 500px; height: auto;"}

#### Riemann solver
---

## Temporal Discretization

#### Explicit time stepping
---

#### _IMEX_ Implicit-Explicit time stepping

#### Implicit time stepping
---

## Newton-Krylov algorithm

