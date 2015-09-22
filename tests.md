---
layout: page
title: Tests
---
$$\mathtt{grim}$$ has been tested extensively in both linear and nonlinear as
well as in special relativistic, and general relativistic regimes.

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

#### Orzag-Tang
