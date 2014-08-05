#include "physics.h"

void setGamma(const struct geometry geom,
              struct fluidElement elem)
{
    elem.gamma = 
        sqrt(1 + geom.gCov[1][1]*elem.primVars[U1]*elem.primVars[U1]
               + geom.gCov[2][2]*elem.primVars[U2]*elem.primVars[U2] 
               + geom.gCov[3][3]*elem.primVars[U3]*elem.primVars[U3]

             + 2*(  geom.gCov[1][2]*elem.primVars[U1]*elem.primVars[U2]
                  + geom.gCov[1][3]*elem.primVars[U1]*elem.primVars[U3] 
                  + geom.gCov[2][3]*elem.primVars[U2]*elem.primVars[U3]
                 )
            );
}

void setDGammaDt(const struct geometry geom,
                 const struct fluidElement elem,
                 struct dFluidElementDt dElemDt)
{
    dElemDt.dGammaDt = 
        ((  geom.gCov[1][1]*elem.primVars[U1]*dElemDt.dPrimVarsDt[U1]
          + geom.gCov[2][2]*elem.primVars[U2]*dElemDt.dPrimVarsDt[U2]
          + geom.gCov[3][3]*elem.primVars[U3]*dElemDt.dPrimVarsDt[U3]
         )
          + (  geom.gCov[1][2]*dElemDt.dPrimVarsDt[U1]*elem.primVars[U2]
             + geom.gCov[1][2]*elem.primVars[U1]*dElemDt.dPrimVarsDt[U2]
             
             + geom.gCov[1][3]*dElemDt.dPrimVarsDt[U1]*elem.primVars[U3]
               geom.gCov[1][3]*elem.primVars[U1]*dElemDt.dPrimVarsDt[U3]

             + geom.gCov[2][3]*dElemDt.dPrimVarsDt[U2]*elem.primVars[U3]
               geom.gCov[2][3]*elem.primVars[U2]*dElemDt.dPrimVarsDt[U3]
            )
        )/elem.gamma;
}

void setUCon(const struct geometry geom,
             struct fluidElement elem)
{
    elem.uCon[0] = elem.gamma/geom.alpha;

    for (int i=1; i<NDIM; i++)
    {
        elem.uCon[i] =   elem.primVars[UU+i] 
                       - elem.gamma*geom.gCon[0][i]*geom.alpha;
    }
}

void setDUConDt(const struct geometry geom,
                const struct fluidElement elem,
                struct dFluidElementDt dElemDt)
{
    dElemDt.dUConDt[0] = dElemDt.dGammaDt/geom.alpha;

    for (int i=1; i<NDIM; i++)
    {
        dElemDt.dUConDt[i] =   dElemDt.dPrimVarsDt[UU+i]
                             - dElemDt.dGammaDt*geom.gCon[0][i]*geom.alpha;
    }

}

void setBCon(const struct geometry geom,
             struct fluidElement elem)
{
    REAL uCov[NDIM];

    conToCov(elem.uCon, geom, uCov);

    elem.bCon[0] =   elem.primVars[B1]*elem.uCov[1]
                   + elem.primVars[B2]*elem.uCov[2] 
                   + elem.primVars[B3]*elem.uCov[3];
    
    for (int i=1; i<NDIM; i++)
    {
        elem.bCon[i] = (  elem.primVars[U3+i] 
                        + elem.bCon[0]*elem.uCon[i]
                       )/elem.uCon[0];
    }
}

void setDBConDt(const struct geometry geom,
                const struct fluidElement elem,
                struct dFluidElementDt dElemDt)
{
    REAL uCov[NDIM], dUCovDt[NDIM];

    conToCov(elem.uCon, geom, uCov);
    conToCov(dElemDt.dUConDt, geom, dUCovDt);

    dElemDt.dBConDt[0] =   elem.primVars[B1]*dUCovDt[1]
                         + dElemDt.dPrimVarsDt[B1]*uCov[1];

                         + elem.primVars[B2]*dUCovDt[2]
                         + dElemDt.dPrimVarsDt[B2]*uCov[2]

                         + elem.primVars[B3]*dUCovDt[3]
                         + dElemDt.dPrimVarsDt[B3]*uCov[3];

    for (int i=1; i<NDIM; i++)
    {
        dElemDt.dBConDt[i] = -(  elem.primVars[U3+i] 
                               + elem.bCon[0]*elem.uCon[i]
                              )*dElemDt.dUConDt[0]/\
                              (elem.uCon[0]*elem.uCon[0])

                             +(  dElemDt.dPrimVarsDt[U3+i] 
                               + dElemDt.dBConDt[0]*elem.uCon[i] 
                               + elem.bCon[0]*dElemDt.dUConDt[i]
                              )/elem.uCon[0];
    }
}

void setFluidElement(const REAL primVars[DOF],
                     const struct geometry geom,
                     struct fluidElement elem)
{
    for (int var=0; var<DOF; var++)
    {
        elem.primVars[var] = primVars[var];
    }

    /* Need to be set in exactly the following order because of the dependecy
       structure */
    setGamma(geom, elem);
    setUCon(geom, elem);
    setBCon(geom, elem);
}

void setDFluidElementDt(const struct fluidElement,
                        const REAL dPrimVarsDt[DOF],
                        const struct geometry geom,
                        struct dFluidElementDt dElemDt)
{
    for (int var=0; var<DOF; var++)
    {
        dElemDt.dPrimVarsDt[var] = dPrimVarsDt[var];
    }

    setDGammaDt(geom, elem, dElemDt);
    setDUConDt(geom, elem, dElemDt);
    setDBConDt(geom, elem, dElemDt);
}
                     
void matterCurrent(const struct fluidElement elem,
                   const struct geometry geom,
                   REAL NUp[NDIM])
{
    for (int mu=0; mu<NDIM; mu++)
    {
        NUp[mu] = elem.primVars[RHO]*elem.uCon[mu];
    }
}

void stressTensor(const struct fluidElement elem,
                  const struct geometry geom,
                  REAL TUpDown[NDIM][NDIM])
{
    REAL pressure = (ADIABATIC_INDEX - 1.)*elem.primVars[UU];
    REAL bCov[NDIM], uCov[NDIM], bSqr;
    
    conToCov(elem.bCon, geom, bCov);
    bSqr = covDotCon(bCov, elem.bCon);

    conToCov(elem.uCon, geom, uCov);

#define DELTA(mu, nu) (mu==nu ? 1 : 0)

    for (int mu=0; mu<NDIM; mu++)
    {
        for (int nu=0; nu<NDIM; nu++)
        {
            TUpDown[mu][nu] =   (  elem.primVars[RHO] + elem.primVars[UU] 
                                 + pressure + bSqr
                                )*elem.uCon[mu]*uCov[nu]

                              + (pressure + 0.5*bSqr)*DELTA(mu, nu)

                              - elem.bCon[mu]*bCov[nu];
        }
    }

#undef DELTA
}


void computeFluxesAndConservedVars(const struct fluidElement elem,
                                   const struct geometry geom,
                                   const int dir,
                                   REAL fluxes[DOF],
                                   REAL conservedVars[DOF])
{
    REAL NUp[NDIM], TUpDown[NDIM][NDIM];
    REAL g = sqrt(-geom.gDet);

    matterCurrent(elem, geom, NUp);
    stressTensor(elem, geom, TUpDown);

    
    fluxes[RHO] = g*NUp[dir];

    fluxes[UU] = g*TUpDown[dir][0];
    fluxes[U1] = g*TUpDown[dir][1];
    fluxes[U2] = g*TUpDown[dir][2];
    fluxes[U3] = g*TUpDown[dir][3];

    fluxes[B1] = g*(elem.bCon[1]*elem.uCon[dir] - elem.bCon[dir]*elem.uCon[1]);
    fluxes[B2] = g*(elem.bCon[2]*elem.uCon[dir] - elem.bCon[dir]*elem.uCon[2]);
    fluxes[B3] = g*(elem.bCon[3]*elem.uCon[dir] - elem.bCon[dir]*elem.uCon[3]);

    conservedVars[RHO] = g*NUp[0];

    conservedVars[UU] = g*TUpDown[0][0];
    conservedVars[U1] = g*TUpDown[0][1];
    conservedVars[U2] = g*TUpDown[0][2];
    conservedVars[U3] = g*TUpDown[0][3];

    conservedVars[B1] = g*(  elem.bCon[1]*elem.uCon[0] 
                           - elem.bCon[0]*elem.uCon[1]);

    conservedVars[B2] = g*(  elem.bCon[2]*elem.uCon[0] 
                           - elem.bCon[0]*elem.uCon[2]);

    conservedVars[B3] = g*(  elem.bCon[3]*elem.uCon[0] 
                           - elem.bCon[0]*elem.uCon[3]);
}

void computeDUDt(const struct fluidElement elem,
                 const struct dFluidElementDt dElemDt,
                 const struct geometry geom,
                 REAL dUDt[DOF])
{
    REAL pressure, dPressureDt, bSqr, dBSqrDt;
    REAL g, tmp1, dTmp1Dt, tmp2, dTmp2Dt;

    REAL bCov[NDIM];
    conToCov(elem.bCon, geom, bCov);
    bSqr = conDotCov(elem.bCon, bCov);

    REAL dBCovDt[NDIM];
    conToCov(dElemDt.dBConDt, geom, dBCovDt);
    
    dBSqrDt = 0.;
    for (int mu=0; mu<NDIM; mu++)
    {
        dBSqrDt += elem.bCon[mu]*dBCovDt[mu] + dElemDt.dBConDt[mu]*bCov[mu];
    }

    pressure = (ADIABATIC_INDEX-1.)*elem.primVars[UU];
    dPressureDt = (ADIABATIC_INDEX-1.)*dElemDt.dPrimVarsDt[UU];
    
    tmp1 = pressure + elem.primVars[RHO] + elem.primVars[UU] + bSqr;
    dTmp1Dt =   dPressureDt + dElemDt.dPrimVarsDt[RHO] 
              + dElemDt.dPrimVarsDt[UU] + dBSqrDt;

    tmp2 = pressure + 0.5*bSqr;
    dTmp2Dt = dPressureDt + 0.5*dBSqrDt;
    
    g = sqrt(geom.gDet);
    
    dUDt[RHO] = g*(  dElemDt.dPrimVarsDt[RHO]*elem.uCon[0]
                   + elem.primVars[RHO]*dElemDt.dUConDt[0]
                  );

#define DELTA(mu, nu) (mu==nu ? 1 : 0)
    for (int var=UU; var<=U3; var++)
    {
        dUDt[var] = g*(  dTmp1Dt*elem.uCon[0]*uCov[var-UU]
                       + tmp1*(  dElemDt.dUConDt[0]*uCov[var-UU]
                               + elem.uCon[0]*dUCovDt[var-UU]
                              ) + dTmp2Dt*DELTA(var, UU)
                       - (  dElemDt.dBConDt[0]*bCov[var-UU]
                          + elem.bCon[0]*dBCovDt[var-UU]
                         )
                      );
    }
#undef DELTA

    dUDt[B1] = g*dElemDt.dPrimVarsDt[B1];
    dUDt[B2] = g*dElemDt.dPrimVarsDt[B2];
    dUDt[B3] = g*dElemDt.dPrimVarsDt[B3];
}

