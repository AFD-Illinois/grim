#ifndef GRIM_INPUT_H_
#define GRIM_INPUT_H_

#define COMPUTE_DIM         (2)
#define NDIM                (4)

/* Domain inputs */
#define N1                  (64)
#define N2                  (64)

/* Number of ghost zones */
#define NG                  (3)

/* Choose the metric here:
 * KERRSCHILD or MINKOWSKI
*/
#define METRIC (MINKOWSKI)
#define M_PI (3.14159265359)
#define H_SLOPE (0.3)
#define A_SPIN (0.9)

#define EPS (1e-5) /* Constant needed for numerical differentiation of metric
                    * to compute the connection coefficients */

/* Boundary conditions: Can be OUTFLOW, MIRROR or PERIODIC */
#define PHYSICAL_BOUNDARY_LEFT_EDGE   (PERIODIC)
#define PHYSICAL_BOUNDARY_RIGHT_EDGE  (PERIODIC)
#define PHYSICAL_BOUNDARY_TOP_EDGE   (PERIODIC)
#define PHYSICAL_BOUNDARY_BOTTOM_EDGE  (PERIODIC)

#define X1_A (0.)
#define X2_A (0.)
#define X1_B (0.)
#define X2_B (0.)

#define TIME_STEPPING EXPLICIT

#define ADIABATIC_INDEX (4./3)
#define CONDUCTION      (0)

/* Double or float precision? */
#define REAL  double

/* Use as ARRAY_ARGS ARRAY_SIZE while declaring arrays in functions arguments */
/* Must use -std=c99 while compiling */
#define ARRAY_ARGS const restrict static

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


#endif /* GRIM_INPUT_H_ */
