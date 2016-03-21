#ifndef GRIM_TORUS_HPP
#define GRIM_TORUS_HPP

#include "../problem.hpp"

/* Internal functions */
double lFishboneMoncrief(double a, double r, double theta);

double lnOfhTerm1(double a,
                double r, double theta, 
                double l);

double lnOfhTerm2(double a,
                double r, double theta, 
                double l);

double lnOfhTerm3(double a,
                double r, double theta, 
                double l);

double computeLnOfh(double a, double r, double theta);

double computeDelta(double a, double r, double theta);

double computeSigma(double a, double r, double theta);

double computeA(double a, double r, double theta);

void applyFloor(grid* prim, fluidElement* elem, geometry* geom,grid* XCoords, int &numReads,int &numWrites);

#endif
