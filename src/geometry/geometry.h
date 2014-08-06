#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include "inputs.h"

/* struct geometry
 *
 * Data structure containing all the geometry related data for a particular
 * zone. All quantities specified in X^mu coordinates. We could be more general
 * and specify all the quantities in x^mu coordinates and numerically transform
 * them to X^mu coordinates during the evolution but its an unnecessary expense.
 * We simply note that whenever xToX() and XTox() has been changed, the
 * function gCov specifying the geometry also should be changed appropriately:
 *
 * ds^2 = g_mu_nu dx^mu dx^nu
 *      = (g_mu_nu (dx^mu/dX^alpha) (dx^nu/dX^beta)) dX^alpha dX^beta
 *      = g_alpha_beta dX^alpha dX^beta
 *
 * Note that we do not store the connection coefficients because they are only
 * used to compute the source terms at the center of a zone and this is only
 * done once inside the zone.
 *
 * The values of this struct are all set by calling setGeometry().
*/
struct geometry
{
    REAL alpha, gDet;
    REAL gCov[NDIM][NDIM];
    REAL gCon[NDIM][NDIM];
};


/* User functions */
void setGeometry(const REAL X[NDIM], 
                 struct geometry* restrict geom);

void covToCon(const REAL vecCov[NDIM],
              const struct geometry* restrict geom,
              REAL vecCon[NDIM]);

void conToCov(const REAL vecCon[NDIM],
              const struct geometry* restrict geom,
              REAL vecCov[NDIM]);

REAL covDotCon(const REAL vecCov[NDIM],
               const REAL vecCon[NDIM]);

REAL gammaDownDownDown(const int eta,
                       const int mu, 
                       const int nu,
                       const REAL X[NDIM]);

/* Functions called by internally inside the geometry module */
void gCovFunc(const REAL X[NDIM], REAL gCov[NDIM][NDIM]);

void gDetFunc(const REAL gCov[NDIM][NDIM], REAL gDet);

void gConFunc(const REAL gCov[NDIM][NDIM],
              const REAL gDet,
              REAL gCon[NDIM][NDIM]);

void XTox(const REAL X[NDIM], REAL x[NDIM]);

#endif /* GRIM_GEOMETRY_H_ */
