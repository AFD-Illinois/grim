#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../geometry/geometry.hpp"
#include "../reconstruction/reconstruction.hpp"

inline int DELTA(int const &mu, int const &nu)
{
  return (mu==nu ? 1 : 0);
}

class fluidElement
{
  void computeEMHDGradients(const double dX[3],
                            int &numReads,
                            int &numWrites
                           );     
  public:
    array zero, one;

    /* fluidElement parameters */
    array tau, chi_emhd, nu_emhd;
    
    array rho, u, u1, u2, u3, B1, B2, B3;
    array pressure, temperature;
    array qTilde, deltaPTilde;
    array q, deltaP;
  
    array gammaLorentzFactor, uCon[NDIM], uCov[NDIM];
    array bSqr, bCon[NDIM], bCov[NDIM];
    array soundSpeed;
    
    array NUp[NDIM];
    array TUpDown[NDIM][NDIM];

    array gradT[NDIM];
    array dtuCov[NDIM];
    array graduCov[NDIM][NDIM];
    array divuCov;
    array deltaP0;
    array q0;
    array bNorm;

    fluidElement(const grid &prim,
                 geometry &geom,
                 int &numReads,
                 int &numWrites
                );
    ~fluidElement();

    void set(const grid &prim,
             geometry &geom,
             int &numReads,
             int &numWrites
            );
    void setFluidElementParameters();
    void computeFluxes(const int direction,
                       grid &flux,
                       int &numReads,
                       int &numWrites
                      );                                

    void computeMinMaxCharSpeeds(const int dir,
                                 array &MinSpeed,
                                 array &MaxSpeed,
                                 int &numReads,
                                 int &numWrites
                                );

    void computeTimeDerivSources(const fluidElement &elemOld,
                                 const fluidElement &elemNew,
                                 const double dt,
                                 grid &sources,
                                 int &numReads,
                                 int &numWrites
                                );

    void computeImplicitSources(grid &sources,
                                array &tauDamp,
                                int &numReads,
                                int &numWrites
                               );

    void computeExplicitSources(const double dX[3],
                                grid &sources,
                                int &numReads,
                                int &numWrites
                               );
    
    void constructTetrads();
    void tetradConToCoordCon(const array vTetrad[NDIM], array vCoord[NDIM]);
    void coordConToTetradCon(const array vCoord[NDIM], array vTetrad[NDIM]);
                             
  private:
    array eCon[NDIM][NDIM];
    array eCov[NDIM][NDIM];
    geometry *geom;
    
    void normalize(array vCon[NDIM]
                  );
    void projectOut(const array vConB[NDIM], 
                    array vConA[NDIM]
                   );
    void normalizeNull(array vCon[NDIM]
                      );
};

class riemannSolver
{
  public:
    fluidElement *elemFace;

    grid *fluxLeft, *fluxRight;
    grid *consLeft, *consRight;

    array minSpeedLeft,  maxSpeedLeft;
    array minSpeedRight, maxSpeedRight;

    riemannSolver(const grid &prim, geometry &geom);
    ~riemannSolver();

    void solve(const grid &primLeft,
               const grid &primRight,
               geometry &geomLeft,
               geometry &geomRight,
               const int dir,
               grid &flux,
               int &numReads,
               int &numWrites
              );
};

#endif /* GRIM_PHYSICS_H_ */
