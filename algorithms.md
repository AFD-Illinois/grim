---
layout: page
title: Algorithms
---

$$\mathtt{grim}$$ uses a second order finite volume method to solve partial
differential equations in the following conservative form
\begin{align}
\frac{\partial U}{\partial t} + \nabla\cdot F & = S
\end{align}
in a spatial volume $$d$$ with the boundary $$\partial d$$, where the
_conserved_ variables $$U \equiv U(P)$$, the _fluxes_ $$F \equiv F(P)$$, and the
_sources_ $$S \equiv S(P)$$ are all functions of the primitive variables to be
solved for, $$P$$. 

#### PDE evolution as a _root finding_ problem
---
In $$\mathtt{grim}$$, the time evolution of the above PDE is cast as a massive
nonlinear root finding problem, by re-writing it as<br>
$$\begin{align}
\mathbf{R}(P^{n+1}) = \frac{\partial U}{\partial t} + \nabla\cdot F - S
\end{align}$$<br>
where $$P^{n+1}$$ are the primitive variables that need to be solved for, and
$$\mathbf{R}(P^{n+1})$$ is a vector of the $$residuals$$ at every spatial
location in the discrete domain. The time evolution of the PDE then is a
question of what the _roots_ $$P^{n+1}$$ of the above set of equations are. The
exact form of $$\partial_t U \equiv (\partial_t U)(P^{n+1}, P^n)$$, $$
\nabla\cdot F \equiv (\nabla \cdot F)(P^{n+1}, P^{n})$$ and $$S \equiv
S(P^{n+1}, P^{n})$$ depend on the chosen temporal discretization scheme,
and are described in the temporal discretization section. The roots are obtained
by the _Newton-Krylov_ algorithm, which only requires the residuals $$R$$ as
inputs. The Jacobian of the system needed for the nonlinear root finding is
assembled automatically with the input residuals.

The residuals $$\mathbf{R}$$ are assembled by discretizing the above
conservative equations using the finite volume formulation, where all the
individual steps are detailed in later sections.

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
where $$\bar{U} = (\int U \Delta v)/\int \Delta v$$, $$\bar{S} = (\int S \Delta
v)/\int \Delta v$$, $$\bar{F}^1 = (\int F^1 dX^2)/\int dX^2$$, and $$\bar{F}^2 =
(\int F^2 dX^1)/\int dX^1$$. The locations $$right$$, $$left$$, $$top$$ and $$bottom$$ are
illustrated below for an individual grid zone. In $$\mathtt{grim}$$, we denote
centers of the grid zones by half indices and faces by integer indices as
shown.<br>

![gridzone](../gridzone.png){: style="max-width: 500px; height: auto;"}

#### Temporal discretization
---
Multiplying further by $$dt$$, and performing the temporal integral $$\int dt
\partial_t \bar{U}$$ over a discrete time interval $$\Delta t$$ gives
$$\begin{align}
\bar{U}_{n+1} - \bar{U}_n + \frac{\int dt \bar{F}^1_{right} - \int
dt\bar{F}^1_{left}}{\Delta X^1} + \frac{\int dt\bar{F}^2_{top} - \int
dt\bar{F}^2_{bottom}}{\Delta X^2} = \int dt \bar{S}
\end{align}$$<br>
where the indices $$n$$, and $$n+1$$ indicate the discrete time levels. The
volume integrals in $$\bar{U}$$, $$\bar{F}$$ and the surface integrals in
$$\bar{F}^{1,2}$$ are evaluated using a second order numerical quadrature, and
thus $$\int dX^1 (.) \rightarrow \Delta X^1 (.)_{i+1/2}$$, and $$\int dX^2 (.)
\rightarrow \Delta X^2 (.)_{j+1/2}$$. For the temporal integrals $$\int dt
\bar{F}^{1,2}$$ and $$\int dt \bar{S}$$, $$\mathtt{grim}$$ can use any of the
following three schemes:<br>

#### (1) Explicit time stepping:

   
  * Fluxes : $$\int dt\bar{F}^{1,2} \rightarrow \Delta t\bar{F}^{1,2}_{n+1/2}$$
  * Sources : $$\int dt\bar{S} \rightarrow \Delta t\bar{S}_{n+1/2}$$
  
  The temporal half index $$n+1/2$$ indicates a half time step. This leads to
  the following discrete equations<br>
   $$\begin{align} \frac{U_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} - U_{n, i+\frac{1}{2}, j+\frac{1}{2}}}{\Delta t} & \\ \nonumber + \frac{F^1_{n+\frac{1}{2}, i+1, j+\frac{1}{2}} - F^1_{n+\frac{1}{2}, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j+1} - F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j}}{\Delta     X^2} & = S_{n+\frac{1}{2}, i+\frac{1}{2}, j+\frac{1}{2}} \end{align}$$

   A complete time step for this scheme is performed in two stages. A _half
   step_ $$n \rightarrow n+1/2$$, to solve for the primitive variables at the
   half step $$P_{n+1/2}$$, which are then used to compute
   $$\bar{F}^{1,2}_{n+1/2} \equiv \bar{F}^{1,2}(P_{n+1/2})$$ and
   $$\bar{S}_{n+1/2} \equiv \bar{S}(P_{n+1/2})$$.  These are then used to perform a
   full step $$n\rightarrow n+1$$, thus completing the time integration over
   $$\Delta t$$.
   
#### (2) _IMEX_ Implicit-Explicit time stepping

  * Fluxes (explicit) : $$\int dt\bar{F}^{1,2} \rightarrow \Delta t\bar{F}^{1,2}_{n+1/2}$$
  * Sources (implicit) : $$\int dt\bar{S} \rightarrow 0.5\Delta t(\bar{S}_{n+1} +
    \bar{S}_n)$$

  This scheme is designed to handle the presence of stiff sources, by treating
  them implicitly, while treating the flux terms explicitly. It leads to the
  following discrete equations<br>
   $$\begin{align} \frac{U_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} - U_{n, i+\frac{1}{2}, j+\frac{1}{2}}}{\Delta t} & \\ \nonumber + \frac{F^1_{n+\frac{1}{2}, i+1, j+\frac{1}{2}} - F^1_{n+\frac{1}{2}, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j+1} - F^2_{n+\frac{1}{2}, i+\frac{1}{2}, j}}{\Delta     X^2} & = 0.5\left(S_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} + S_{n, i+\frac{1}{2}, j+\frac{1}{2}}\right) \end{align}$$

   The computation of the fluxes at the half step is the same as in the explicit
   scheme.
   
#### (3) Implicit time stepping:
  
  * Fluxes : $$\int dt\bar{F}^{1,2} \rightarrow 0.5\Delta t(\bar{F}^{1,2}_{n+1}  + \bar{F}^{1,2}_n)$$.
  * Sources : $$\int dt\bar{S} \rightarrow 0.5\Delta t(\bar{S}_{n+1} +
    \bar{S}_n)$$

    This scheme is much more expensive than the explicit and imex schemes as
    will be described in the solver section.  However, a fully implicit scheme has
    no Courant limits and is useful for testing new physics when the
    characteristics are not known accurately. The discrete equations are<br>
    $$\begin{align} \label{eq:fvm_implicit_time_stepping}\frac{U_{n+1,
    i+\frac{1}{2}, j+\frac{1}{2}} - U_{n, i+\frac{1}{2}, j+\frac{1}{2}}}{\Delta
    t} & \\ \nonumber+ 0.5\left(\frac{F^1_{n+1, i+1, j+\frac{1}{2}} - F^1_{n+1, i,j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n+1, i+\frac{1}{2}, j+1} - F^2_{n+1,i+\frac{1}{2}, j}}{\Delta X^2}\right) \\ \nonumber + 0.5\left(\frac{F^1_{n, i+1,j+\frac{1}{2}} - F^1_{n, i, j+\frac{1}{2}}}{\Delta X^1} + \frac{F^2_{n,i+\frac{1}{2}, j+1} - F^2_{n, i+\frac{1}{2}, j}}{\Delta X^2}\right) & = 0.5\left(S_{n+1, i+\frac{1}{2}, j+\frac{1}{2}} + S_{n, i+\frac{1}{2}, j+\frac{1}{2}}\right)\end{align}$$


#### Spatio-temporal derivatives in the source terms
---
The volume averaged source terms may contain spatio-temporal derivatives
$$\bar{S} \equiv \bar{S}(\partial_t, \partial_i)$$. Whenever present, the
temporal derivatives in the source terms are approximated as $$\partial_t(.)
\approx ((.)_{n+1} - (.)_n)/\Delta t$$, and the spatial derivatives are
approximated using _slope limited_ derivatives. The derivative of a quantity in
the $$i$$ zone is computed from the values at the neighbouring points using
$$\partial_{x}(.)\approx limiter\left[(.)_{i-1}, (.)_i, (.)_{i+1}\right] +
error$$, where $$error \sim O(\Delta x^2)$$ in smooth flows and $$O(\Delta x)$$
in the presence of discontinuities. The spatial derivatives of the required
quantities are computed using the primitive variables at the $$n+1/2$$ half-step
for the explicit and imex schemes, whereas it is a combination of the $$n$$ and
$$n+1$$ steps for the implicit scheme.

#### Reconstruction
---
The finite volume method evolves the zone-averaged conserved variables
$$\bar{U}_n \rightarrow \bar{U}_{n+1}$$, where $$\bar{U}_{n,n+1} \approx U_{n,
n+1, i+1/2, j+1/2} \equiv U(P_{n, n+1, i+1/2, j+1/2})$$. Therefore, the
conserved variables $$U$$ and the primitive variables $$P$$ are located at zone
centers $$(i+1/2, j+1/2)$$, whereas the fluxes $$F(P)$$ need to be computed at
the $$right$$, $$left$$, $$top$$ and $$bottom$$ face centers . Thus, the need to
_reconstruct_ the values of the primitive variables $$P_{i+1/2, j+1/2}$$ to the
zone faces. This is performed using a reconstruction operator $$R$$ which takes
in primitive variables within a certain radius and constructs a polynomial
interpolant from which the edge states can be computed. To achieve an error of
$$O(\Delta x^2)$$, a linear interpolant is sufficient, for which the
reconstruction operator has a stencil width of three grid zones. The operator
can act in two directions depending on input order, $$R^+_{i+1/2} \equiv
R(P_{i-1/2}, P_{i+1/2}, P_{i+3/2})$$ whose output is the $$right$$ state
$$P_{i+1}$$, and $$R^-_{i+1/2} \equiv R(P_{i+3/2}, P_{i+1/2}, P_{i-1/2})$$,
whose output is the $$left$$ state $$P_{i}$$.

Illustrated below are the sequence of steps needed to compute $$F^1_{left}$$ and
$$F^1_{right}$$.

![reconstruction](../reconstruction.png){: style="max-width: 550px; height: auto;"}

  * (a) Apply the operator $$R^+_{i-1/2} \equiv R(P_{i-3/2}, P_{i-1/2},
      P_{i+1/2})$$ to compute the primitive variables $$P^-_{i}$$ at the left
      side of the $$left$$ face (index = $$i$$).

  * (b) $$R^-_{i+1/2} \equiv R(P_{i+3/2}, P_{i+1/2},
      P_{i-1/2})$$ to compute the primitive variables $$P^+_{i}$$ at the right
      side of the $$left$$ face (index = $$i$$).

    $$R^+_{i+1/2} \equiv R(P_{i-1/2}, P_{i+1/2}, P_{i+3/2})$$ to compute the
    primitive variables $$P^-_{i+1}$$ at the left side of the $$right$$ face
    (index = $$i+1$$).

  * (c) $$R^-_{i+3/2} \equiv R(P_{i+5/2}, P_{i+3/2},
      P_{i+1/2})$$ to compute the primitive variables $$P^+_{i+1}$$ at the right
      side of the $$right$$ face (index = $$i+1$$).

After the reconstruction procedure described in the above sequence of steps, we
have the primitive variables at the left $$P^-_i$$ and the right $$P^+_i$$ sides
of the $$left$$ face and at the left $$P^-_{i+1}$$ and the right $$P^+_{i+1}$$
sides of the $$right$$ face. The fluxes at each face are then a function of the
primitive variables at either side of the face, $$F^1_{left} \equiv
F^1(P^-_i, P^+_i)$$ and $$F^1_{right} \equiv F^1(P^-_{i+1}, P^+_{i+1})$$. These
are then computed using the _Riemann solver_.

The above procedure outlines the reconstruction in one-dimension. Since
$$\mathtt{grim}$$ uses a structured mesh, the same one-dimension reconstruction
operator along with the accompanying sequence of steps are followed to compute
the edge states for the $$top$$ and $$bottom$$ faces.

#### Riemann solver
---

Given the left and right states on either side of a face, an approximate Riemann
solver is used to compute the flux at the face. Shown below are the primitive
variables after reconstruction.

![riemann_problem](../riemann_problem.png){: style="max-width: 300px; height: auto;"}

$$\mathtt{grim}$$ uses the Lax-Friedrichs flux for its approximate Riemann
solver, which requires as an input the maximum characteristic speed $$c_{max}$$
of the physical model being solved. The Lax-Friedrichs flux is given by,<br>
$$\begin{align}
F^1_i = \frac{1}{2} (F^1(P^+_i) + F^1(P^-_i)) - c_{max}(U(P^+_i) - U(P^-_i))
\end{align}$$
