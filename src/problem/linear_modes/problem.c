#include "../problem.h"

#if (CONDUCTION)
  REAL kappaProblem;
  REAL tauProblem; 

  void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1])
  {
    elem->kappa = kappaProblem;
    elem->tau   = tauProblem;

  #if (MODE==FULL_EMHD_2D)
  REAL Rho = elem->primVars[RHO];
//  if(Rho<RHO_FLOOR_MIN)
//     Rho=RHO_FLOOR_MIN;
  REAL U = elem->primVars[UU];
//  if(U<UU_FLOOR_MIN)
//     U = UU_FLOOR_MIN;
//
  REAL P   = (ADIABATIC_INDEX-1.)*U;
  REAL T   = P/Rho;
//  if(T<1.e-12)
//    T=1.e-12;
  
  REAL cs  = sqrt(  ADIABATIC_INDEX*P
                  / (Rho + (ADIABATIC_INDEX*U))
                 );

  REAL beta = 1./(Rho* 1 * cs*cs*T);
    
  REAL tau    = 1.;
  elem->kappa = tau/beta/T;
  //elem->kappa = 0.1;
  elem->tau   = tau; 

  #endif
  }
#endif

#if (VISCOSITY)
  REAL etaProblem;
  REAL tauVisProblem; 

  void setViscosityParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1])
  {
    elem->eta     = etaProblem;
    elem->tauVis  = tauVisProblem;

  #if (MODE==FULL_EMHD_2D)
    REAL Rho = elem->primVars[RHO];
//    if(Rho<RHO_FLOOR_MIN)
//      Rho=RHO_FLOOR_MIN;
    REAL U   = elem->primVars[UU];
//    if(U<UU_FLOOR_MIN)
//      U = UU_FLOOR_MIN;
    REAL P   = (ADIABATIC_INDEX-1.)*U;
    REAL T   = P/Rho;
//    if(T<1.e-12)
//      T=1.e-12;
    REAL cs  = sqrt(  ADIABATIC_INDEX*P
                  / (Rho + (ADIABATIC_INDEX*U))
                 );
    
    REAL beta = 0.5/(1 * cs*cs*Rho);

    REAL tau     = 1.;
    elem->eta    = 0.5*tau/beta;
    //elem->eta    = 0.1;
    elem->tauVis = tau;

  #endif
  }
#endif

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  REAL randNum;
  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);


  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);

      #if (MODE==FULL_EMHD_2D) 
      /* Eigenvalue = -0.5533585207638141 - 3.6262571286888425*I
       *
       * kappa = rho * 1 * cs**2 * tau
       * eta   = rho * 1 * cs**2 * tau
       * tau   = 1
       * Gamma = 4/3
       * k1 = 2*pi
       * k2 = 4*pi
       *
       * For this problem, the values of kappa and eta are set directly in
       * setConductionParameters() and setViscosityParameters()
       *
       *('Eigenvalue   = ', -0.5533585207638108 - 3.626257128688849*I)

        (delta_rho, ' = ', -0.5185225240822464 - 0.1792647678001874*I)
        (delta_u, ' = ', 0.551617073639382)
        (delta_u1, ' = ', 0.008463122479547853 + 0.011862022608466373*I)
        (delta_u2, ' = ', -0.16175466371870748 - 0.034828080823603495*I)
        (delta_u3, ' = ', 1.1150888301102058e-17 + 2.3900040672376563e-17*I)
        (delta_B1, ' = ', -0.059737949796407556 - 0.0335170750615094*I)
        (delta_B2, ' = ', 0.029868974898203757 + 0.01675853753075467*I)
        (delta_B3, ' = ', -1.317599308756542e-17 - 3.254129323467699e-17*I)
        (delta_q, ' = ', 0.5233486841539429 + 0.04767672501939605*I)
        (delta_dp, ' = ', -0.2909106062057659 - 0.021594520553365606*I)
       *
       * */
        
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.1;
        primVars0[B2]  = 0.3;
        primVars0[B3]  = 0.;
        primVars0[PHI] = 0.;
        primVars0[PSI] = 0.;

        deltaPrimVars[RHO] = -0.518522524082246 - 0.1792647678001878*I;
        deltaPrimVars[UU]  = 0.5516170736393813;
        deltaPrimVars[U1]  = 0.008463122479547856 + 0.011862022608466367*I;
        deltaPrimVars[U2]  = -0.16175466371870734 - 0.034828080823603294*I;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = -0.05973794979640743 - 0.03351707506150924*I;
        deltaPrimVars[B2]  = 0.02986897489820372 + 0.016758537530754618*I;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = 0.5233486841539436 + 0.04767672501939603*I;
        /* IMPORTANT NOTE: Balbusaur outputs deltaP, but psi = -deltaP */
        deltaPrimVars[PSI] = 0.2909106062057657 + 0.02159452055336572*I;

//        deltaPrimVars[RHO] =  0.8245277155993238;
//        deltaPrimVars[UU]  = -0.2655263062409864 + 0.11295573948161661*I;
//        deltaPrimVars[U1]  =  0.01870577997997143 - 0.07463791562430622*I;
//        deltaPrimVars[U2]  =  0.13548478237175182 + 0.0029582968057185223*I;
//        deltaPrimVars[U3]  =  0.;
//        deltaPrimVars[B1]  =  0.013767510382234251 + 0.13241920582661107*I;
//        deltaPrimVars[B2]  = -0.006883755191117101 - 0.06620960291330555*I;
//        deltaPrimVars[B3]  =  0.;
//        deltaPrimVars[PHI] = -0.38685433587246587 + 0.11071197829033866*I;
//        deltaPrimVars[PSI] =  -0.1599729569151908 + 0.054267588035646866*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 4*M_PI;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }


      #elif (MODE==HYDRO_SOUND_MODE_1D) /* Eigenvalue = 3.09362659024*I */
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;

        deltaPrimVars[RHO] = 0.345991032308;
        deltaPrimVars[UU]  = 0.922642752822;
        deltaPrimVars[U1]  = -0.170354208129;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }

      #elif (MODE==ENTROPY_WAVE_1D) /* Eigenvalue = 0 */
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;

        deltaPrimVars[RHO] = 1.;
        deltaPrimVars[UU]  = 0.;
        deltaPrimVars[U1]  = 0.;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }

      #elif (MODE==CONDUCTION_STABLE_1D) 
        /* Eigenvalue = -0.498597173331 - 0.857832357798*I
         * kappa = 0.1
         * tau = 1.01818181818182 */
        kappaProblem = .1;
        tauProblem   = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.0000000001;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
        primVars0[PHI] = 0.;
        
        deltaPrimVars[RHO] = 0.911376933183;
        deltaPrimVars[UU]  = 0.030751595371 - 0.0635975709194*I;
        deltaPrimVars[U1]  = 0.124428706971 - 0.0723215917578*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = -0.332658158111 + 0.181734443922*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==CONDUCTION_STABLE_2D) 
        /* Eigenvalue = -0.498689052044 - 1.23434614343*I
         * kappa = 0.1
         * tau = 1.01818181818182 */
        kappaProblem = .1;
        tauProblem   = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.02;
        primVars0[B3]  = 0.;
        primVars0[PHI] = 0.;

        deltaPrimVars[RHO] = 0.912960868047;
        deltaPrimVars[UU]  = 0.0441633411305 - 0.0470501442451*I;
        deltaPrimVars[U1]  = 0.068161459988 - 0.0280266780212*I;
        deltaPrimVars[U2]  = 0.111191793414 - 0.0444339558103*I;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = -0.00130516937615 + 6.41592876069e-05*I;
        deltaPrimVars[B2]  = 0.00130516937615 - 6.41592876069e-05*I;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = -0.352802085664 + 0.134521891027*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 2*M_PI;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==VISCOSITY_2D) 
        /* Eigenvalue = -0.488186438695 + 0.516141145988*I
         * eta = 0.1
         * tau = 1.01818181818182 */
        etaProblem    = .1;
        tauVisProblem = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
        primVars0[PSI] = 0.;

        deltaPrimVars[RHO] = 0.34457248379 - 8.32667268469e-17*I;
        deltaPrimVars[UU]  = 0.918859956773;
        deltaPrimVars[U1]  = -0.180848290246 - 0.00323785674963*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PSI] = -0.0624512526648 - 0.0186932205343*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==VISCOSITY_1D)
        etaProblem    = .1;
        tauVisProblem = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.02;
        primVars0[B3]  = 0.;
        primVars0[PSI] = 0.;

        deltaPrimVars[RHO] = 0.34591909261 + 8.32667268469e-17*I;
        deltaPrimVars[UU]  = 0.922450913626;
        deltaPrimVars[U1]  = -0.170824276817 - 0.00015896242418*I;
        deltaPrimVars[U2]  = 0.00288321673818 + 0.000950778108505*I;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.00697678484543 + 1.91989053629e-05*I;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PSI] = 0.0129157109207 + 0.00431579622043*I;

        REAL k1 = 2*M_PI;

        REAL complex mode = cexp(I*k1*XCoords[1]);

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==ALFVEN_2D) 
        /* Eigenvalue = -0.488186438695 + 0.516141145988*I
         * eta = 0.1
         * tau = 1.01818181818182 */
        
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
	
        deltaPrimVars[RHO] = 0.;
        deltaPrimVars[UU]  = 0.;
        deltaPrimVars[U1]  = 0.;
        deltaPrimVars[U2]  = 1.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 1.00005;
        deltaPrimVars[B3]  = 0.;

	etaProblem    = .1;
        tauVisProblem = 1.01818181818182;
	primVars0[PSI]  = 0.;
	deltaPrimVars[PSI]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==FIREHOSE) 
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = .001;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
	primVars0[PSI]  = 0.;
	
        deltaPrimVars[RHO] = 0.345907557128 + 2.08166817117e-17*I;
        deltaPrimVars[UU]  = 0.922420152341;
        deltaPrimVars[U1]  = -0.171584412179 - 4.06388721081e-05*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
	deltaPrimVars[PSI]  = -0.00691108899638 - 0.000221744119852*I;

	etaProblem    = .1;
        tauVisProblem = 10.;

	REAL FireSeed = 1.e-8;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
	  if(var!=U2 && var!=B2)
	    INDEX_PETSC(primOldGlobal, &zone, var) =  
	      primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
	mode = cexp(I*(k1*XCoords[1]*20 + k2*XCoords[2]) );
	INDEX_PETSC(primOldGlobal, &zone, U2) = creal(FireSeed*mode);
	INDEX_PETSC(primOldGlobal, &zone, B2) = -creal(FireSeed*mode*sqrt(11./3));
	//PetscRandomGetValue(randNumGen, &randNum);
	//INDEX_PETSC(primOldGlobal, &zone, U2) = FireSeed*randNum;
	//PetscRandomGetValue(randNumGen, &randNum);
        //INDEX_PETSC(primOldGlobal, &zone, B2) = FireSeed*randNum;
      #endif

    }
  }
  
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);
}


void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR)
    {
        primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR)
    {
        primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR;
    }
  }
}

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
                                     const int X1Start, const int X2Start,
                                     const int X1Size, const int X2Size,
                                     REAL primTile[ARRAY_ARGS TILE_SIZE])
{

}

void applyAdditionalProblemSpecificBCs
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL primTile[ARRAY_ARGS TILE_SIZE]
)
{

}

void applyProblemSpecificFluxFilter
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
)
{

}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}


void writeProblemSpecificData(PetscViewer parametersViewer,
  const struct problemData problemSpecificData[ARRAY_ARGS 1]) {

}
