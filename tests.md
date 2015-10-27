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
built on top of the $$\mathtt{sagemath}$$ computer algebra system, which takes
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
vector $$k_x = 2\pi, k_y = 4 \pi$$, which is misaligned with the background
magnetic field $$B_{x0} = 0.1, B_{y0} = 0.3$$. Evidently, both the wave
propagation vector $$\mathbf{k}$$ and the background magnetic field
$$\mathbf{B}_0$$ are misaligned with the numerical grid.

The initial conditions are set with the following eigenmode

| Variable     | Background state | Perturbed value                                  |
|:------------:|:----------------:|:--------------------------------------------------------------------:|
|$$\rho$$      |1.                |-0.518522524082246 - 0.1792647678001878*$$i$$     | 
|$$u$$         |1.                |0.5516170736393813| 
|$$u^1$$       |0.                |0.008463122479547856 + 0.011862022608466367*$$i$$ | 
|$$u^2$$       |0.                |-0.16175466371870734 - 0.034828080823603294*$$i$$ | 
|$$u^3$$       |0.                |0.                                                | 
|$$B^1$$       |0.1               | -0.05973794979640743 - 0.03351707506150924*$$i$$ | 
|$$B^2$$       |0.3               |0.02986897489820372 + 0.016758537530754618*$$i$$  | 
|$$B^3$$       |0.                |0.                                                | 
|$$q$$         |0.                |0.5233486841539429 + 0.04767672501939605*$$i$$    | 
|$$\Delta P$$  |0.                |0.2909106062057659 + 0.021594520553365606*$$i$$   | 

All the variables are 

To run this test  
  * Open `src/CMakeLists.txt`  
  * Set the problem with `set(PROBLEM "linear_modes")`  
  * Open `src/problem/linear_modes/CMakeLists.txt`  
  * Set this mode using `set(PROBLEM "FULL_EMHD_2D")`  

The plot below shows that all evolved variables in $$\mathtt{grim}$$ converge to their
respective analytic solutions at the expected second order.
![linear_modes_convergence](../linear_modes_convergence_full_emhd.png){:style="max-width: 500px; height: auto;"}

[^balbusaur]: [$$\mathtt{balbusaur}$$](https://cloud.sagemath.com/projects/4aac3c0b-be00-496a-a5c9-9b6079c8ab02/files/grim/src/problem/linear_modes/Balbusaur.sagews): A framework for automated linear analysis. Hosted on [SageMathCloud](http://www.sagemath.com)

#### EMHD Shock solutions

![rho_viscosity](../stationary_shock_rho_viscosity.png){:style="max-width:600px; height: auto;"}
![dP_viscosity](../stationary_shock_psi_viscosity.png){:style="max-width:600px; height: auto;"}

#### Hydrostatic Conducting Atmosphere

![atmosphere_heat_flux](../atmosphere_test_heat_flux.png){:style="max-width:600px; height: auto;"}
![convergence_atmosphere](../convergence_test_atmosphere.png){:style="max-width:600px; height: auto;"}

#### Bondi Inflow without backreaction

![bondi_heat_flux](../bondi_heat_flux.png){:style="max-width:600px; height: auto;"}
![bondi_pressure_anisotropy](../bondi_inflow_viscosity_soln.png){:style="max-width:550px; height: auto;"}

#### Anisotropic Conduction _Snake_ Test

<iframe width="550" height="550" frameborder="0" src="https://app.box.com/embed/preview/wzz038mcelmflq8xvu8g?view=&sort=&direction=ASC&theme=light"></iframe>

#### Magneto Thermal Instability

<iframe width="550" height="550" frameborder="0" src="https://app.box.com/embed/preview/65l28xjkqoqi1yw1r29h?view=&sort=&direction=ASC&theme=light"></iframe>

![radial_mag_field_growth](../MTI_growth_radial_magnetic_energy.png){:style="max-width:600px; height: auto;"}

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

