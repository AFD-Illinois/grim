---
layout: page
title: Algorithms
---

$$\mathtt{grim}$$ uses a second order finite volume method to solve equations in the following conservative form
\begin{align}
\frac{\partial U}{\partial t} + \nabla\cdot F & = S
\end{align}
in a spatial volume $$d$$ with the boundary $$\partial d$$.

## Spatial Discretization
The spatial volume $$d$$ is discretized by transforming into a set of chosen computational coordinates $$X^\mu$$, in which the coordinate axes are aligned to the transformed boundary $$\partial D$$. The domain in these coordinates is then discretized using a structured uniform hexahedral mesh.

#### Grid generation
---
![grid](../grid.png){: style="max-width: 500px; height: auto;"}

#### Grid indexing
---
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

