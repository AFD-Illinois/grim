---
layout: page
title: Tests
---
$$\mathtt{grim}$$ has been tested extensively in linear, nonlinear, special
and general relativistic regimes. The tests below are grouped according to the
physical model being solved, with subsections describing individual tests.

## Extended MHD

#### Linear modes
It is important to check that any numerical implementation of the EMHD model can
reproduce the corresponding linear theory with an error that falls off at the
expected order of spatio-temporal discretization. In order to perform this
test, one needs to first derive the linear theory of the EMHD model.

The governing equations of the EMHD model are considerably more complex than
those of ideal MHD. In particular, the inclusion of both anisotropic pressure
and conduction, which are sourced by spatio-temporal derivatives projected along
the magnetic field lines make it challenging to derive the linear theory of this
model. This is especially so for arbitrary inclinations of the wave vectors
$$\mathbf{k}$$ with the magnetic field $$\mathbf{B}$$. Even in special cases
where $$\mathbf{k}\parallel\mathbf{B}$$ and $$\mathbf{k} \perp \mathbf{B}$$, the
derivation is prone to errors if done manually. To address this issue, we have
written a general linear analysis package $$\mathtt{balbusaur}$$[^balbusaur],
built on top of the $$\mathtt{sagemath}$$[^sagemath] computer algebra system, which takes
as input the governing equations of any model, and generates the characteristic
matrix of the corresponding linear theory. The eigenvectors of this matrix are
then used as initial conditions in $$\mathtt{grim}$$, and their numerical
evolutions checked against the corresponding analytic solutions.

With $$\mathtt{balbusaur}$$, we are able to compute the linear modes of the
EMHD model for a relativistically hot configuration where rest mass energy
density  $$\rho$$ is comparable to the internal energy density $$u$$, including
both anisotropic pressure and conduction, and for any inclinations of the wave
vector $$\mathbf{k}$$ with the magnetic field $$\mathbf{B}$$.  

To thoroughly test the numerical implementation, we choose a mode which excites
the variables $$\{\rho, u, u^1, u^2, B^1, B^2, q, \Delta P\}$$, with a wave
vector $$k_x = 2\pi, k_y = 4 \pi$$, that is misaligned with the background
magnetic field $$B_{x0} = 0.1, B_{y0} = 0.3$$. Evidently, both the wave
propagation vector $$\mathbf{k}$$ and the background magnetic field
$$\mathbf{B}_0$$ are misaligned with the numerical grid.

The initial conditions are set with the following eigenmode having eigenvalue $$\omega$$ = -0.5533585207638108 - 3.626257128688849$$i$$.

| Variable     | Background state | Perturbed value                                  |
|:------------:|:----------------:|:------------------------------------------------:|
|$$\rho$$      |1.                |-0.518522524082246-0.1792647678001878$$i$$     | 
|$$u$$         |1.                |0.5516170736393813                                | 
|$$u^1$$       |0.                |0.008463122479547856+0.011862022608466367$$i$$ | 
|$$u^2$$       |0.                |-0.16175466371870734-0.034828080823603294$$i$$ | 
|$$u^3$$       |0.                |0.                                                | 
|$$B^1$$       |0.1               | -0.05973794979640743-0.03351707506150924$$i$$ | 
|$$B^2$$       |0.3               |0.02986897489820372+0.016758537530754618$$i$$  | 
|$$B^3$$       |0.                |0.                                                | 
|$$q$$         |0.                |0.5233486841539429+0.04767672501939605$$i$$    | 
|$$\Delta P$$  |0.                |0.2909106062057659+0.021594520553365606$$i$$   | 

All the variables are initialized as $$x = x_0 + A\delta x$$, where $$x$$ is the variable, $$x_0$$ is the background state, $$\delta x$$ is the perturbation and $$A$$ is the amplitude.

To run this test  
  * Open `src/CMakeLists.txt`  
  * Set the problem with `set(PROBLEM "linear_modes")`  
  * Open `src/problem/linear_modes/CMakeLists.txt`  
  * Set this mode using `set(PROBLEM "FULL_EMHD_2D")`  

The test has been run with the following options in `src/problem/linear_modes/CMakeLists.txt`

    # Time stepping options: EXPLICIT, IMEX or IMPLICIT
    set(TIME_STEPPING "IMEX")
    set(DT "0.01")
    set(DT_DUMP ".1")
    set(START_TIME "0.")
    set(FINAL_TIME ".5")
    set(START_DUMP_COUNTER "0")
    set(COURANT ".9")
    set(MAX_DT_INCREMENT "1.3")

    # Domain size. If the problem is 1D then N2 is ignored.
    set(N1 "512")
    set(N2 "512")

    # Geometry
    set(METRIC "MINKOWSKI")
    set(EPS    "1e-5")

    # Domain
    set(X1_A  "0.")
    set(X1_B  "1.")
    set(X2_A  "0.")
    set(X2_B  "1.")

    # Boundary conditions
    set(PHYSICAL_BOUNDARY_LEFT_EDGE   "PERIODIC")
    set(PHYSICAL_BOUNDARY_RIGHT_EDGE  "PERIODIC")
    set(PHYSICAL_BOUNDARY_TOP_EDGE    "PERIODIC")
    set(PHYSICAL_BOUNDARY_BOTTOM_EDGE "PERIODIC")

    # Reconstrution options
    # MP5, MONOTONIZED_CENTRAL or MIN_MOD
    set(RECONSTRUCTION "MONOTONIZED_CENTRAL")

    # Initial condition parameters 
    # Mode options: 
    # 1) ENTROPY_WAVE_1D
    # 2) HYDRO_SOUND_MODE_1D
    # 3) CONDUCTION_STABLE_1D
    # 4) CONDUCTION_STABLE_2D
    # 5) VISCOSITY_2D
    # 6) VISCOSITY_1D
    # 7) ALFVEN_2D
    # 8) FIREHOSE
    # 9) FULL_EMHD_2D
    set(AMPLITUDE "1e-8")
    set(MODE      "FULL_EMHD_2D")

The plot below shows that all evolved variables in $$\mathtt{grim}$$ converge to their
respective analytic solutions at the expected second order.
![linear_modes_convergence](../linear_modes_convergence_full_emhd.png){:style="max-width: 500px; height: auto;"}

[^balbusaur]: [$$\mathtt{balbusaur}$$](https://cloud.sagemath.com/projects/4aac3c0b-be00-496a-a5c9-9b6079c8ab02/files/grim/src/problem/linear_modes/Balbusaur.sagews): A framework for automated linear analysis. Hosted on [SageMathCloud](http://www.sagemath.com)
[^sagemath]: Sage Mathematics Software Version 6.7

#### EMHD Shock solutions

The presence of sufficient viscosity can smoothen a shock and connect the left
and the right sides with a well-defined solution. The hyperbolic nature of the
dissipation in the EMHD model leads to new features in the shock structure
which have been qualitatively described in Chandra et. al. Here, we solve the
shock structure in the EMHD model as a boundary value problem with the left and
the right states fixed to their values given by the Rankine-Hugoniot jump
conditions. We then use this as a reference solution to check the EMHD shock
solutions obtained from $$\mathtt{grim}$$, which solves the EMHD equations as an
initial value problem. 

The boundary value solutions are obtained using a global Newton root finder. We
are looking for a steady state nonlinear solution of the EMHD equations and
hence set the time derivatives $$\partial_t \rightarrow 0$$. Since we are
interested in the $$\mathit{continuous}$$ shock sub-structure, we approximate all
spatial derivatives $$\partial_x$$ by central differences with a truncation error
$$O(\Delta x^8)$$. Thus we have a set of coupled discrete nonlinear equations
$$R(P_i) = 0$$, where $$P_i$$ are the primitive variables at $$i=0, 1, ..., N_x$$ and
$$N_x$$ is the chosen spatial resolution of the numerical grid. The system is
iterated upon starting from a smooth initial guess using the Newton's method
combined with a numerical Jacobian assembled to machine precision. The
iterations are continued till we achieve machine precision error $$O(10^{-14})$$.

We find that the major contribution to the shock structure comes from the
pressure anisotropy with the role of the heat conduction in determining the
values of the thermodynamic variables inside the shock being marginal. The EMHD
equations have higher order corrections $$\sim q u^\mu \nabla_mu(\tau_R/(\chi
P^2)), \Delta P u^\mu \nabla_\mu(\tau_R/(\rho \nu P))$$ that we expect to
contribute in strong nonlinear regimes, and indeed we see that the shock
structure differs as we turn on and turn off these terms.

There exists an upper limit to the strength of the shock that can be solved for
using the EMHD model. Higher mach number shocks require a larger viscosity (or
$$\Delta P$$) to smoothly connect the left and the right states. Beyond a certain
critical value of $$\Delta P$$, the theory loses hyperbolicity, and eventually
causality and stability. The root of this problem lies in the fact that
ultimately, the theory is a second order perturbation $$\sim q^2, \Delta P^2$$,
about an equilibrium and as the dissipative effects become stronger, the
validity of the expansion breaks down. We obtain the thresholds of this
breakdown by expanding about an non-equilibrium using a background heat flux
$$q_0$$ and assembling the characteristic matrix using $$\mathtt{balbusaur}$$. As $$q_0$$
is increased, hyperbolicity is first broken upon the appearance of imaginary
values in the eigenvalues. As $$q_0$$ is further increased, the real values of
the eigenvalues exceed $$c$$, indicating a breakdown of causality. Finally, there
is a loss of stability with the linear modes growing exponentially with the
fastest growth at $$k \rightarrow \infty$$.

![rho_viscosity](../bvp_solver_vs_grim_with_conduction_256_chi_5.png){:style="max-width:600px; height: auto;"}
![rho_viscosity](../bvp_solver_HO_vs_No_HO.png){:style="max-width:600px; height: auto;"}
![dP_viscosity](../bvp_solver_residuals_with_conduction_256_chi_5_random_initial_guess.png){:style="max-width:600px; height: auto;"}

#### Hydrostatic Conducting Atmosphere
Relativisitic conduction in curved space-time contains qualitatively new
features when compared to non-relativistic conduction because the heat flux is
driven by _red-shifted_ temperature gradients $$\nabla_\mu T + T a_\mu$$ where
$$a_\mu = u^\nu \nabla_\mu u_\nu$$ is the four-acceleration. We test this effect
with a hydrostatic fluid configuration in a Schwarzschild metric.

![atmosphere_heat_flux](../atmosphere_test_heat_flux.png){:style="max-width:600px; height: auto;"}
![convergence_atmosphere](../convergence_test_atmosphere.png){:style="max-width:600px; height: auto;"}

#### Bondi Inflow without backreaction

![bondi_heat_flux](../bondi_heat_flux.png){:style="max-width:600px; height: auto;"}
![bondi_pressure_anisotropy](../bondi_inflow_viscosity_soln.png){:style="max-width:550px; height: auto;"}

#### Anisotropic Conduction _Snake_ Test

<iframe width="550" height="550" frameborder="0" src="https://app.box.com/embed/preview/wzz038mcelmflq8xvu8g?view=&sort=&direction=ASC&theme=light"></iframe>

#### Magneto Thermal Instability

<iframe width="600" height="550" frameborder="0" src="https://app.box.com/embed/preview/65l28xjkqoqi1yw1r29h?view=&sort=&direction=ASC&theme=light"></iframe>

![radial_mag_field_growth](../MTI_growth_radial_magnetic_energy.png){:style="max-width:600px; height: auto;"}

#### Heat Flux Driven buoyancy Instability

<iframe width="560" height="525" frameborder="0" src="https://app.box.com/embed/preview/f7ce5l1dzbj3o7nmgi6ov1rqirrz6t9k?direction=ASC&theme=light"></iframe>

![theta_mag_field_growth](../hbi_growth.png){:style="max-width:600px; height: auto;"}

#### Firehose Instability

<iframe width="560" height="550" frameborder="0" src="https://app.box.com/embed/preview/uuk907qk75cz0imd34szod4d8uvequ48?direction=ASC&theme=light"></iframe>

![firehose_growth](../firehose_unstable_mode_growth.png){:style="max-width:600px; height: auto;"}

## Ideal MHD tests

#### Sound modes

#### Komissarov shock tests
 Relativistic MHD shock test suite from Komissarov (1999)[^komissarov_1999] and
 reproduced with corrections by Gammie, McKinney and Toth[^harm_paper]. The
 tests involve starting with _left_ and _right_ states for all thermodynamic
 variables and letting the system evolve to a desired time. The solution is then
 compared to the known solutions in the literature. The initial states for each
 test are tabulated below. To run this test suite

  * Open `src/CMakeLists.txt`
  * Set the problem with `set(PROBLEM "shock_tests")`
  * Open `src/problem/shock_tests/CMakeLists.txt`
  * Choose the specific shock test to run, for example `set(PROBLEM "FAST_SHOCK")`

 The test have all been run with the following options in `src/problem/shock_tests/CMakeLists.txt`

    # Time stepping options: EXPLICIT, IMEX or IMPLICIT
    set(TIME_STEPPING "IMEX")
    set(DT "1e-5")
    set(START_TIME "0.")
    set(FINAL_TIME "100.")
    set(COURANT "0.5")
    set(MAX_DT_INCREMENT "1.3")
    
    # Domain size. If the problem is 1D then N2 is ignored.
    set(COMPUTE_DIM "1")
    set(N1 "512")
    
    # Physics variables
    set(ADIABATIC_INDEX "4./3")
    set(CONDUCTION "OFF")
    set(VISCOSITY  "ON")
    
    # Geometry
    set(METRIC "MINKOWSKI")
    set(EPS    "1e-5")
    
    # Domain
    set(X1_A  "-4.")
    set(X1_B  "4.")
    
    # Boundary conditions
    set(PHYSICAL_BOUNDARY_LEFT_EDGE   "OUTFLOW")
    set(PHYSICAL_BOUNDARY_RIGHT_EDGE  "OUTFLOW")
    
    # Reconstrution options
    # MP5, MONOTONIZED_CENTRAL or MIN_MOD
    set(RECONSTRUCTION "MONOTONIZED_CENTRAL")

[^komissarov_1999]: _A Godunov-type scheme for relativistic magnetohydrodynamics_, Komissarov (1999)
[^harm_paper]: _HARM: A Numerical Scheme for General Relativistic Magnetohydrodynamics_, Gammie, McKinney, and Toth (2003)

###### Fast shock

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.          |25.48        | 
|$$Pressure$$  |1.          |367.5        | 
|$$u^1$$       |25.         |1.091        | 
|$$u^2$$       |0.          |0.3923       | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |20.         |20.          | 
|$$B^2$$       |25.02       |49.          | 
|$$B^3$$       |0.          |0.           | 

![fast_shock_rho](../Komissarov_fast_shock_rho.png){: style="max-width: 600px; height: auto;"}
![fast_shock_ux](../Komissarov_fast_shock_ux.png){: style="max-width: 600px; height: auto;"}

###### Slow shock

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.          |3.323        | 
|$$Pressure$$  |10.         |55.36        | 
|$$u^1$$       |1.53        |0.9571       | 
|$$u^2$$       |0.          |-0.6822      | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |10.         |10.          | 
|$$B^2$$       |18.28       |14.49        | 
|$$B^3$$       |0.          |0.           | 

![slow_shock_rho](../Komissarov_slow_shock_rho.png){: style="max-width: 600px; height: auto;"}
![slow_shock_ux](../Komissarov_slow_shock_ux.png){: style="max-width: 600px; height: auto;"}

###### Switch on slow

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.78e-3     |0.01         | 
|$$Pressure$$  |0.1         |1.           | 
|$$u^1$$       |-0.765      |0.           | 
|$$u^2$$       |-1.386      |0.           | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |1.          |1.           | 
|$$B^2$$       |1.022       |0.           | 
|$$B^3$$       |0.          |0.           | 

![switch_on_rho](../Komissarov_switch_on_rho.png){: style="max-width: 600px; height: auto;"}
![switch_on_ux](../Komissarov_switch_on_ux.png){: style="max-width: 600px; height: auto;"}

###### Switch off fast

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |0.1         |0.562        | 
|$$Pressure$$  |1.          |10.          | 
|$$u^1$$       |-2.         |-0.212       | 
|$$u^2$$       |0.          |-0.590       | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |2.          |2.           | 
|$$B^2$$       |0.          |4.71         | 
|$$B^3$$       |0.          |0.           | 

![switch_off_rho](../Komissarov_switch_off_rho.png){: style="max-width: 600px; height: auto;"}
![switch_off_ux](../Komissarov_switch_off_ux.png){: style="max-width: 600px; height: auto;"}

###### Alfven wave

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.          |1.           | 
|$$Pressure$$  |1.          |1.           | 
|$$u^1$$       |0.          |3.7          | 
|$$u^2$$       |0.          |5.76         | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |3.          |3.           | 
|$$B^2$$       |3.          |-6.857       | 
|$$B^3$$       |0.          |0.           | 

###### Shock tube 1

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.          |0.1          | 
|$$Pressure$$  |1000.       |1.           | 
|$$u^1$$       |0.          |0.           | 
|$$u^2$$       |0.          |0.           | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |1.          |1.           | 
|$$B^2$$       |0.          |0.           | 
|$$B^3$$       |0.          |0.           | 

![shock_tube_1_rho](../Komissarov_shock_tube_1_rho.png){: style="max-width: 600px; height: auto;"}
![shock_tube_1_ux](../Komissarov_shock_tube_1_ux.png){: style="max-width: 600px; height: auto;"}

###### Shock tube 2

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.          |0.1          | 
|$$Pressure$$  |30.         |1.           | 
|$$u^1$$       |0.          |0.           | 
|$$u^2$$       |0.          |0.           | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |0.          |0.           | 
|$$B^2$$       |20.         |0.           | 
|$$B^3$$       |0.          |0.           | 

![shock_tube_2_rho](../Komissarov_shock_tube_2_rho.png){: style="max-width: 600px; height: auto;"}
![shock_tube_2_ux](../Komissarov_shock_tube_2_ux.png){: style="max-width: 600px; height: auto;"}


###### Collision

| Variable     | Left state | Right state |
|:------------:|:----------:|:-----------:|
|$$\rho$$      |1.          |1.           | 
|$$Pressure$$  |1.          |1.           | 
|$$u^1$$       |5.          |-5.          | 
|$$u^2$$       |0.          |0.           | 
|$$u^3$$       |0.          |0.           | 
|$$B^1$$       |10.         |10.          | 
|$$B^2$$       |10.         |-10.         | 
|$$B^3$$       |0.          |0.           | 

![collision_rho](../Komissarov_collision_rho.png){: style="max-width: 600px; height: auto;"}
![collision_ux](../Komissarov_collision_ux.png){: style="max-width: 600px; height: auto;"}

#### Orzag-Tang

<iframe width="550" height="550" frameborder="0" src="https://app.box.com/embed/preview/h3bwc7i0iib9s4sm3ljzb2h3nhucwuk2?view=&sort=&direction=ASC&theme=light"></iframe>

#### Komissarov Magnetized Cylindrical Explosion

<iframe width="550" height="550" frameborder="0" src="https://app.box.com/embed/preview/yonoxdr9wrrw48cm88tc4owy5vrgf9xt?view=&sort=&direction=ASC&theme=light"></iframe>

