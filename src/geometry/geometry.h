#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include <math.h>
#include "../inputs.h"
#include "macros.h"

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
 * The values of this struct are all set by calling setGeometry().
*/
struct geometry
{
  REAL XCoords[NDIM];
  REAL alpha, gDet;
  REAL gCov[NDIM][NDIM];
  REAL gCon[NDIM][NDIM];
};


/* User functions */
void setGeometry(const REAL X[ARRAY_ARGS NDIM], 
                 struct geometry geom[ARRAY_ARGS 1]);

void covToCon(const REAL vecCov[ARRAY_ARGS NDIM],
              const struct geometry geom[ARRAY_ARGS 1],
              REAL vecCon[ARRAY_ARGS NDIM]);

void conToCov(const REAL vecCon[ARRAY_ARGS NDIM],
              const struct geometry geom[ARRAY_ARGS 1],
              REAL vecCov[ARRAY_ARGS NDIM]);

REAL covDotCon(const REAL vecCov[ARRAY_ARGS NDIM],
               const REAL vecCon[ARRAY_ARGS NDIM]);

REAL gammaDownDownDown(const int eta,
                       const int mu, 
                       const int nu,
                       const REAL X[ARRAY_ARGS NDIM]);

void XTox(const REAL X[ARRAY_ARGS NDIM], 
          REAL x[ARRAY_ARGS NDIM]);

/* Functions called by internally inside the geometry module */
void gCovFunc(const REAL X[ARRAY_ARGS NDIM],
              REAL gCov[ARRAY_ARGS NDIM][NDIM]);

REAL gDetFunc(REAL gCov[ARRAY_ARGS NDIM][NDIM]);

void gConFunc(REAL gCov[ARRAY_ARGS NDIM][NDIM],
              const REAL gDet,
              REAL gCon[ARRAY_ARGS NDIM][NDIM]);

#endif /* GRIM_GEOMETRY_H_ */
