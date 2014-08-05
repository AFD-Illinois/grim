#ifndef GRIM_INPUT_H_
#define GRIM_INPUT_H_

#define COMPUTE_DIM         (2)
#define NDIM                (4)
#define N1                  (64)
#define N2                  (64)
#define NG                  (2)

/* Choose the geometry here:
 * KERRSCHILD or MINKOWSKI
*/
#define KERRSCHILD

/* Double or float precision? */
#define REAL                (double)

/* The entire global domain is divided into tiles that are small enough to fit
 * into the cache of the compute node or an accelerator. This technique
 * optimizes cache usage and makes explicit use of the deep memory hierarchies
 * which are prevalent in all kinds of processors and accelerators. 
 * 
 * Caution: 1) The tile sizes need to divide the global domain size in each
 *             direction exactly!
 *          2) We want a tile size that is as large as the cache. Too large a
 *             size and the code can fail silently because local memory in the
 *             OpenCL specification is a different amount for different 
 *             machines.
 *             NEED TO PUT A CHECK FOR THIS.
 *
*/
#define TILE_SIZE_X1        (8)
#define TILE_SIZE_X2        (8)

#ifdef KERRSCHILD
#define A_SPIN              (0.9375)    // Black hole spin
#define ADIABATIC_INDEX     (4./3.)

// Floor values:
#define SMALL               (1e-15)     // Minimum value of rho and u allowed when
                                        // computing fluxes inside the residual
                                        // evaluation function. Set to the limits of
                                        // REAL precision.

#define RHO_FLOOR_0         (1e-4)      // rho_floor = RHO_FLOOR_0 * radius^RHO_FLOOR_EXPONENT
#define RHO_FLOOR_EXPONENT  (-1.5)           
#define UU_FLOOR_0          (1e-6)      // u_floor = UU_FLOOR_0 * radius^UU_FLOOR_EXPONENT
#define UU_FLOOR_EXPONENT   (-2.5)      
#define GAMMA_MAX           (50.)       // Code cannot handle arbitrary large lorentz factors. 
                                        // Limit the Lorentz factor to this value.
// End of floor values

// Initial conditions:
// Parameters for Fishbone-Moncrief solution: 
// R_DISK_INNER_EDGE, R_PRESSURE_MAX, ENTROPY_CONSTANT
#define R_DISK_INNER_EDGE   (6.)        // Inner edge of the disk in the initial
                                        // conditions

#define R_PRESSURE_MAX      (12.)       // Pressure max of the disk in the initial
                                        // conditions

#define ENTROPY_CONSTANT    (1e-3)      // Constant entropy factor in polytropic 
                                        // EOS for the initial conditions
                                        // P = ENTROPY_CONSTANT rho^ADIABATIC_INDEX

#define BETA                (1e2)       // Plasma beta for the initial conditions
// End of initial conditions

/**
 * Boundaries for the computational domain X^mu = {t, X1, X2, phi}. The
 * physical domain is set by the transformation in XTox(). All
 * computational operations are performed in X^mu coordinates. The
 * discretization is uniform in X^mu coordinates and the solution for
 * all variables are in X^mu coordinates.
 *
 *   +----> X1
 *   |
 *   |
 *   v
 *   X2
 *      (X1_A, X2_A)     (X1_B, X2_B)
 *           +----------------+
 *           |                |
 *           |                |
 *           |                |
 *           |                | 
 *           +----------------+
 *      (X1_D, X2_D)     (X1_C, X2_C)
 *
*/
#define R_A     (.98*(1. + sqrt(1. - A_SPIN*A_SPIN))) 
#define X1_A    (log(R_A)) /* X1_A as a function of R_A must be consistent with
                              functional dependence set in xtoX() */
#define X2_A    (1e-10)

#define R_B     (40.)
#define X1_B    (log(R_B))
#define X2_B    (1e-10)

#define R_C     (40.)
#define X1_C    (log(R_C))
#define X2_C    (M_PI - 1e-10)

#define R_D     (.98*(1. + sqrt(1. - A_SPIN*A_SPIN)))
#define X1_D    (log(R_D))
#define X2_D    (M_PI - 1e-10)
/* End of boundaries for the computational domain.*/

/** Transformation parameters between x^mu and X^mu: 
 *  1) H_SLOPE: 1 <= H_SLOPE <= 0.
 *     More points near midplane as H_SLOPE -> 0
*/
#define H_SLOPE             (0.3)
/* End of transformation parameters. */


#endif /* KERRSCHILD */

/* Numerical differencing parameter for calculating the connection */
#define EPS (1e-5)

/* Variable mnemonics */
#define RHO (0)
#define UU (1)
#define U1 (2)
#define U2 (3)
#define U3 (4)
#define B1 (5)
#define B2 (6)
#define B3 (7)
#define DOF (8)

/* Boundary mnemonics. See boundary.h for description */
#define TILE_BOUNDARY       (99)
#define OUTFLOW             (100)
#define MIRROR              (101)
#define CONSTANT            (102)
#define PERIODIC            (103)
#define NONE                (104)


// iTile, jTile have ranges [-NG, TILE_SIZE+NG)
#define INDEX_LOCAL(iTile,jTile,var) (iTile+NG + \
                                      (TILE_SIZE_X1+2*NG)*(jTile+NG + \
                                      (TILE_SIZE_X2+2*NG)*(var)))

// i, j have ranges [0, N1), [0, N2)
#define INDEX_GLOBAL(i,j,var) (var + DOF*((i)+(N1)*(j)))

#define i_TO_X1_CENTER(i) (X1_START + (i + 0.5)*DX1)
#define j_TO_X2_CENTER(j) (X2_START + (j + 0.5)*DX2)
#define i_TO_X1_FACE(i) (X1_START + (i)*DX1)
#define j_TO_X2_FACE(j) (X2_START + (j)*DX2)


#endif /* GRIM_INPUT_H_ */
