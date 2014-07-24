#include "gridzone.h"

/**
 *  Convert from physical Kerr-Schild coordinates x^mu = {t, r, theta, phi} to
 *  the computational coordinates X^mu = {t, X1, X2, phi}. This is useful for
 *  concentrating grid points where there are needed. Here we use the following
 *  transformation:
 *
 *  r = exp(X1) - (1)
 *  theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) - (2)
 *
 *  Eqn (1) concentrates points near the horizon and Eqn (2) concentrates points
 *  towards the midplane as the parameter H_SLOPE is decreased from 1 to 0.
 *  
 *  Ref: "A Measurement of the Electromagnetic Luminosity of a Kerr Black Hole"
 *       Jonathan C. McKinney and Charles F. Gammie
 *
 *  @param x Physical coordinates x^mu = {t, r, theta, phi}
 *  @param X Computational coordinates X^mu = {t, X1, X2, phi}
*/
void xToX(const REAL x[NDIM], REAL X[NDIM])
{
    X[0] = x[0];        // t' = t
    X[1] = log(x[1]);   // X1 = log(r)
    X[2] = 0.;          // theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2).
                        // How do I get X2 given a theta?. Set to 0 for now.
    X[3] = x[3];        // phi' = phi
}

/**
 *  Convert from computational coordinates X^mu = {t, X1, X2, phi} to the
 *  physical coordinates x^mu = {t, r, theta, phi}. This is useful for
 *  concentrating grid points where there are needed. Here we use the following
 *  transformation:
 *
 *  r = exp(X1) - (1)
 *  theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) - (2)
 *
 *  Eqn (1) concentrates points near the horizon and Eqn (2) concentrates points
 *  towards the midplane as the parameter H_SLOPE is decreased from 1 to 0.
 *  
 *  Ref: "A Measurement of the Electromagnetic Luminosity of a Kerr Black Hole"
 *       Jonathan C. McKinney and Charles F. Gammie
 *
 *  @param X Computational coordinates X^mu = {t, X1, X2, phi}
 *  @param x Physical coordinates x^mu = {t, r, theta, phi}
*/
void XTox(const REAL X[NDIM], REAL x[NDIM])
{
    x[0] = X[0];        // t = t'
    x[1] = exp(X[1]);   // r = exp(X1)
    x[2] = M_PI*X[2] + 0.5*(1 - H_SLOPE)*sin(2*M_PI*X[2]); 
                        // theta = pi*X2 + 0.5*(1 - H_SLOPE)*X2
    x[3] = X[3];        // phi = phi'
}

/* Get the coordinates corresponding to a location inside a grid zone.
 *
 *         +---> X1
 *         |
 *         |
 *         v
 *         X2
 *         
 *         (i,j)  (i+0.5,j)
 *           +-----F1-----+
 *           |            |
 *           |            |
 * (i,j+0.5) F2    C      |
 *           |            |
 *           |            |
 *           +------------+
 *
 *      C   --  Center (CENTER)
 *      F1  --  Face in the X1 direction (FACE_X1)
 *      F2  --  Face in the X2 direction (FACE_X2)
 *
 * @param input: zone, Zone in which the coordinates are needed
 * @param input: location, Location inside the zone where the coordinates are needed
 * @param output: X, The coordinates X^mu
*/
void getXCoords2D(const struct gridZone2D zone,
                  const int location,
                  REAL X[NDIM])
{
    if (location==CENTER) // (i+0.5,j+0.5)
    {
        X[0] = 0.; // Don't care about the time coordinate. Just set to 0.
        X[1] = X1_A + (zone.i + 0.5)*zone.dX1;
        X[2] = X2_A + (zone.j + 0.5)*zone.dX2;
        X[3] = 0.;
    } 
    else if (location==FACE_X1) // (i+0.5,j)
    {
        X[0] = 0.;
        X[1] = X1_A + (zone.i + 0.5)*zone.dX1;
        X[2] = X2_A + zone.j*zone.dX2;
        X[3] = 0.;
    }
    else if (location==FACE_X2) // (i,j+0.5)
    {
        X[0] = 0.;
        X[1] = X1_A + zone.i*zone.dX1;
        X[2] = X2_A + (zone.j + 0.5)*zone.dX2;
        X[3] = 0.;
    }
    else if (location==CORNER) // (i,j)
    {
        X[0] = 0.;
        X[1] = X1_A + zone.i*zone.dX1;
        X[2] = X2_A + zone.j*zone.dX2;
        X[3] = 0.;
    }
}

int main()
{
    return(0.);
}
