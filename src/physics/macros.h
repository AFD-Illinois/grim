#ifndef GRIM_PHYSICS_MACROS_H_
#define GRIM_PHYSICS_MACROS_H_

/* Primitive variable mnemonics */
#define RHO             (0) /* Rest mass density */
#define UU              (1) /* Internal energy */
#define U1              (2)
#define U2              (3)
#define U3              (4)
#define B1              (5)
#define B2              (6)
#define B3              (7)
#if (CONDUCTION == 1 && VISCOSITY == 1)
  #define PHI           (8)
  #define PSI           (9)
  #define DOF           (10)
#elif (CONDUCTION == 1 && VISCOSITY == 0)
  #define PHI           (8)
  #define DOF           (9)
#elif (CONDUCTION == 0 && VISCOSITY == 1)
  #define PSI           (8)
  #define DOF           (9)
#else
  #define DOF           (8)
#endif
/* DOF == Total Degrees Of Freedom == Total number of equations solved for */

#define DELTA(mu, nu) (mu==nu ? 1 : 0)

/* Indices for elem.moments[] */
#define N_UP(mu) (mu)
#define T_UP_DOWN(mu, nu) (nu + NDIM*(mu) + NDIM)

/* Indices for the Christoffel symbols */
#define GAMMA_UP_DOWN_DOWN(eta,mu,nu) (eta+NDIM*(mu+NDIM*(nu) ) )


#endif
