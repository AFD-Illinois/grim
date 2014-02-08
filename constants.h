#define COMPUTE_DIM 2
#define NDIM 4
#define N1 64
#define N2 64
#define NG 2
#define X1_MIN 0.
#define X1_MAX 1.
#define X2_MIN 0.
#define X2_MAX 1.
#define DX1 (X1_MAX-X1_MIN)/N1
#define DX2 (X2_MAX-X2_MIN)/N2

#define H_SLOPE 0.3
#define R0 0.

#define REAL double
#define TILE_SIZE_X1 16
#define TILE_SIZE_X2 16

#define RHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7
#define DOF 8

// iTile, jTile have ranges [-NG, TILE_SIZE+NG)
#define INDEX_LOCAL(iTile,jTile,var) (iTile+NG + (TILE_SIZE_X1+2*NG)*(jTile+NG + (TILE_SIZE_X2+2*NG)*var))

// i, j have ranges [0, N1), [0, N2)
#define INDEX_GLOBAL(i,j,var) ((i) + N1*((j) + N2*(var)))

