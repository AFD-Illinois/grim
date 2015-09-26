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
The spatial volume $$d$$, is discretized by transforming into a set of chosen
computational coordinates $$X^\mu$$, in which the coordinate axes are aligned to
the transformed boundary $$\partial D$$. The boundaries are then simply defined
by a constant coordinate on all sides of the domain. The domain in these
coordinates is then discretized using a structured uniform hexahedral mesh.
![grid](../grid.png){: style="max-width: 500px; height: auto;"}

A common case involves discretizing spherical domains. This is used for example,
to study accretion flows. To construct a spherical mesh, a computational
coordinate system $$X^\mu \equiv \{A, B, C, D\}$$ is chosen, which relates to the
cartesian $$x^\mu \equiv \{t, x, y, z\}$$ coordinates by<br>
$$\begin{align}
B & = \log(r) \\
\pi C + \frac{1 - h}{2} \sin(2 \pi C) & = \theta
\end{align}$$<br>
where $$r = \sqrt{x^2 + y^2}$$, $$\theta = \tan^{-1}({y/x})$$ and $$C\in [0,
1]$$. With this coordinate mapping, a uniform grid in the $$X^\mu$$ coordinates
leads to an exponential packing ($$dr = e^BdB$$) of the grid zones at the inner
boundary and a concentration of the grid zones at the midplane ($$d\theta =
\pi(1 + (1-h)cos(2\pi C) ) dC$$) which is controlled by the parameter $$h$$.
The discretized domain is illustrated below for $$h = 0.3$$.

![spherical_grid](../spherical_grid_generation.png){: style="max-width: 650px; height: auto;"}

#### Finite Volume Method
---
Muliplying the conservative equations by the volume element of a discrete grid zone
in the $$X^\mu$$ coordinate system, $$\Delta v = dX^1dX^2$$, and using the
divergence theorem leads to,<br>
$$\begin{align}
\partial_t \bar{U} + \frac{\bar{F}^1_{right} - \bar{F}^1_{left}}{\Delta X^1} + \frac{\bar{F}^2_{top} - \bar{F}^2_{bottom}}{\Delta X^2} = \bar{S}
\end{align}$$<br>
where $$\bar{U} = \int U \Delta v/\int \Delta v$$, $$\bar{S} = \int S \Delta v/\int
\Delta v$$ and $$\bar{F}^1 = \int F^1 dX^2/\int dX^2$$, $$\bar{F}^2 = \int F^2
dX^1/\int dX^1$$. The locations $$right$$, $$left$$, $$top$$ and $$bottom$$ are
illustrated below for an individual grid zone. In $$\mathtt{grim}$$, we denote
centers of the grid zones by half indices and faces by integer indices as
shown.<br>

![gridzone](../gridzone.png){: style="max-width: 500px; height: auto;"}

#### Temporal discretization
---
Multiplying by $$dt$$, and performing the
integral over a discrete time interval $$\Delta t$$, $$\int dt \partial_t(.)
\rightarrow (.)_{n+1} - (.)_{n}$$, we get
$$\begin{align}
\bar{U}_{n+1} - \bar{U}_n + \frac{\int dt \bar{F}^1_{right} - \int
dt\bar{F}^1_{left}}{\Delta X^1} + \frac{\int dt\bar{F}^2_{top} - \int
dt\bar{F}^2_{bottom}}{\Delta X^2} = \int dt \bar{S}
\end{align}$$
where the indices $$n$$, and $$n+1$$ indicate the discrete time levels.  We
evaluate the volume and surface integrals using a second order numerical
quadrature, and thus $$\int dX^1 (.) \rightarrow \Delta X^1 (.)_{i+1/2}$$, and
$$\int dX^2 (.) \rightarrow \Delta X^2 (.)_{j+1/2}$$. For the temporal integral
$$\int dt (.)$$, we use three schemes:<br>

#### (1) Explicit time stepping:

 * $$\int dt(.) \rightarrow \Delta t (.)_{n+1/2}$$, where the temporal
   half index $$n+1/2$$ indicates a half time step. This leads to the following
   discrete equations<br>
   $$\begin{align} \frac{U_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} - U_{n, i+\frac{1}{2}, j+\frac{1}{2}}}{\Delta t} & \\ \nonumber + \frac{F^1_{n+\frac{1}{2}, i+1, j+\frac{1}{2}} - F^1_{n+\frac{1}{2}, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j+1} - F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j}}{\Delta     X^2} & = S_{n+\frac{1}{2}, i+\frac{1}{2}, j+\frac{1}{2}} \end{align}$$
   
#### (2) _IMEX_ Implicit-Explicit time stepping

  * $$\int dt(.) \rightarrow \Delta t (.)_{n+1/2}$$ for the integrals involves
    spatial fluxes $$\int dt \bar{F}^1$$, $$\int dt \bar{F}^2$$, as well as for any
    spatial derivatives that may appear in $$\int dt \bar{S}$$. But all the terms
    not involving derivatives in $$\int dt\bar{S}$$ are treated as $$\int dt(.)
    \rightarrow 0.5\Delta t\left((.)_n + (.)_{n+1} \right)$$. Thus we get<br>
   $$\begin{align} \frac{U_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} - U_{n, i+\frac{1}{2}, j+\frac{1}{2}}}{\Delta t} & \\ \nonumber + \frac{F^1_{n+\frac{1}{2}, i+1, j+\frac{1}{2}} - F^1_{n+\frac{1}{2}, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j+1} - F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j}}{\Delta     X^2} & = 0.5\left(S_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} + S_{n, i+\frac{1}{2}, j+\frac{1}{2}}\right) \end{align}$$
   
#### (3) Implicit time stepping:
$$\int dt(.) \rightarrow 0.5\Delta t\left((.)_n + (.)_{n+1} \right)$$. Performing the integrals, with explicit time stepping  we have,
and for the case of implicit time stepping,
$$\begin{align} \label{eq:fvm_implicit_time_stepping}
\frac{U_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} - U_{n, i+\frac{1}{2}, j+\frac{1}{2}}}{\Delta t} & \\ \nonumber
+ 0.5\left(\frac{F^1_{n+1, i+1, j+\frac{1}{2}} - F^1_{n+1, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n+1, i+\frac{1}{2}, j+1} - F^2_{n+1, i+\frac{1}{2}, j}}{\Delta X^2}\right) \\ \nonumber
+ 0.5\left(\frac{F^1_{n, i+1, j+\frac{1}{2}} - F^1_{n, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n, i+\frac{1}{2}, j+1} - F^2_{n, i+\frac{1}{2}, j}}{\Delta X^2}\right)
 & = 0.5\left(S_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} + S_{n, i+\frac{1}{2}, j+\frac{1}{2}}\right)
\end{align}$$

The volume averaged source
terms may contain spatio-temporal derivatives $$\bar{S} \equiv
\bar{S}(\partial_t, \partial_i)$$. All the temporal derivatives 

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

