#ifndef GRIM_RECONSTRUCT_MACROS_H_
#define GRIM_RECONSTRUCT_MACROS_H_

/* RECONSTRUCTION options */
#define MONOTONIZED_CENTRAL (0)
#define MIN_MOD             (1)
#define MP5                 (2)
#define PPM                 (3)

/* The following macros taken from PLUTO. Needed for MP5 */
#define MINMOD(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
#define MY_MAX(a,b)  ( (a) >= (b) ? (a) : (b) ) 
#define MY_MIN(a,b)  ( (a) <= (b) ? (a) : (b) ) 
#define MEDIAN(a,b,c) (a + MINMOD(b-a,c-a))

/* Number of ghost zones */
#if (RECONSTRUCTION==MONOTONIZED_CENTRAL || RECONSTRUCTION==MIN_MOD)
  #define NG  (2)
#elif (RECONSTRUCTION==MP5 || RECONSTRUCTION==PPM)
  #define NG  (3)
#endif

#endif
