#include "grim.h"

#ifdef MINKOWSKI
void SetInitialCondition(TS ts, Vec X)
{
    DM dmda;
    TSGetDM(ts, &dmda);
    REAL ***x;

    DMDAVecGetArrayDOF(dmda, X, &x);

    for (PetscInt j=0; j<N2; j++) {
        for (PetscInt i=0; i<N1; i++) {
            REAL xcoord = DX1/2. + i*DX1;
            REAL ycoord = DX2/2. + j*DX2;

            REAL xcenter = 0.5;
            REAL ycenter = 0.5;

            REAL r = PetscSqrtScalar((xcoord-xcenter)*(xcoord-xcenter)+
                                            (ycoord-ycenter)*(ycoord-ycenter));
            
#ifdef TRANSPORT_TEST
            x[j][i][RHO] = 1. + exp(-r*r/0.01);

            x[j][i][UU] = 1./(ADIABATIC_INDEX-1);
            x[j][i][U1] = 4.95;
            x[j][i][U2] = 4.95;
            x[j][i][U3] = 0.;
            x[j][i][B1] = 0.;
            x[j][i][B2] = 0.;
            x[j][i][B3] = 0.;
#endif

#ifdef ORZAG_TANG_TEST
            x[j][i][RHO] = 25./(36.*M_PI);
            x[j][i][UU] = 5./(12.*M_PI*(ADIABATIC_INDEX - 1.));
            REAL vx = -0.5*PetscSinScalar(2*M_PI*ycoord);
            REAL vy = 0.5*PetscSinScalar(2*M_PI*xcoord);
            REAL vz = 0.;
            REAL gamma = 1./PetscSqrtScalar(1 - vx*vx - vy*vy - vz*vz);
            x[j][i][U1] = gamma*vx;
            x[j][i][U2] = gamma*vy;
            x[j][i][U3] = gamma*vz;
            x[j][i][B1] =
                -1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(2*M_PI*ycoord);
            x[j][i][B2] =
                1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(4*M_PI*xcoord);
            x[j][i][B3] = 0.0;
#endif

#ifdef KOMISSAROV_CYLINDRICAL_EXPLOSION_TEST
            REAL R = 1.;
            REAL rho_outside = 1e-4, rho_inside = 1e-2;
            REAL u_outside = 3e-5/(5./3 - 1); 
            REAL u_inside = 1./(5./3 - 1);
            REAL alpha = 20.;
            REAL norm_factor = (1. - tanh(-R))/2.;

            x[j][i][RHO] = ((rho_inside - rho_outside)*
                            (1. - tanh(pow(r, alpha)-R))
                            /2./norm_factor + rho_outside);

            x[j][i][UU] = ((u_inside - u_outside)*
                           (1. - tanh(pow(r, alpha)-R))
                           /2./norm_factor + u_outside);

            REAL vx = 0.;
            REAL vy = 0.;
            REAL vz = 0.;
            REAL gamma = 1.;
            x[j][i][U1] = gamma*vx;
            x[j][i][U2] = gamma*vy;
            x[j][i][U3] = gamma*vz;
            x[j][i][B1] = .1;
            x[j][i][B2] = 0.;
            x[j][i][B3] = 0.0;
#endif

        }
    }

    DMDAVecRestoreArrayDOF(dmda, X, &x);
}
#endif

void SetInitialCondition(TS ts, Vec Prim)
{
    DM dmda;
    PetscScalar ***prim, r, theta;
    PetscScalar X1, X2;

    PetscScalar l, delta, sigma, A, lnOfh;
    PetscScalar thetaIn, deltaIn, sigmaIn, AIn;
    PetscScalar uPhiTmp, uconBL[NDIM], uconKS[NDIM], uconKSPrime[NDIM]; 
    PetscScalar gcovBL[NDIM][NDIM], gconBL[NDIM][NDIM], transform[NDIM][NDIM];
    PetscScalar rhoMax = 0., uMax = 0., bsqrMax = 0., rhoAv=0.;
    PetscScalar q, norm, betaActual;
    PetscScalar AA, BB, CC, DD, discriminent, mu;

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscScalar randNum;
    PetscRandom randNumGen;
    PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
    PetscRandomSetType(randNumGen, PETSCRAND48);

    TSGetDM(ts, &dmda);

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    PetscScalar AVector[X2Size+2*NG][X1Size+2*NG];

    Vec localPrim;
    DMGetLocalVector(dmda, &localPrim);
    DMDAVecGetArrayDOF(dmda, localPrim, &prim);

	l =  ( ( (pow(A_SPIN, 2) - 2.*A_SPIN*sqrt(R_MAX) + pow(R_MAX, 2.)) *
		     ( (-2.*A_SPIN*R_MAX*(pow(A_SPIN, 2.) - 2.*A_SPIN*sqrt(R_MAX) +
		        pow(R_MAX, 2.))) / sqrt(2.*A_SPIN*sqrt(R_MAX) + (-3. +
                R_MAX)*R_MAX) +
		       ((A_SPIN + (-2. + R_MAX)*sqrt(R_MAX))*(pow(R_MAX, 3) +
				pow(A_SPIN, 2)*(2. + R_MAX))) / sqrt(1 + (2.*A_SPIN) /
			   pow(R_MAX, 1.5) - 3./R_MAX) ) ) / \
           (pow(R_MAX, 3)*sqrt(2.*A_SPIN*sqrt(R_MAX) +
           (-3.+R_MAX)*R_MAX)*(pow(A_SPIN,2) + (-2.+R_MAX)*R_MAX)) );

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++)
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {
            
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            BLCoords(&r, &theta, X1, X2);

            delta = r*r - 2*M*r + A_SPIN*A_SPIN;
            sigma = r*r + A_SPIN*A_SPIN*cos(theta)*cos(theta);
            A = (r*r + A_SPIN*A_SPIN)*(r*r + A_SPIN*A_SPIN) - 
                delta*A_SPIN*A_SPIN*sin(theta)*sin(theta);

            thetaIn = M_PI/2.;

            deltaIn = R_MIN*R_MIN - 2*M*R_MIN + A_SPIN*A_SPIN;
            sigmaIn = R_MIN*R_MIN + A_SPIN*A_SPIN*cos(thetaIn)*cos(thetaIn);
            AIn = (R_MIN*R_MIN + A_SPIN*A_SPIN)*(R_MIN*R_MIN + A_SPIN*A_SPIN) - 
                  deltaIn*A_SPIN*A_SPIN*sin(thetaIn)*sin(thetaIn);

            if (r >=R_MIN) {

			    lnOfh = 0.5*log((1. + sqrt(1. + 4.*(l*l*sigma*sigma)*delta/\
                                (A*sin(theta)*A*sin(theta))))/(sigma*delta/A))-
			            0.5*sqrt(1. + 4.*(l*l*sigma*sigma)*delta /
					             (A*A*sin(theta)*sin(theta))) -2.*A_SPIN*r*l/A-
			           (0.5*log((1. + sqrt(1. +
				                           4.*(l*l*sigmaIn*sigmaIn)*deltaIn/ \
				        (AIn*AIn*sin(thetaIn)*sin(thetaIn)))) /
				        (sigmaIn * deltaIn / AIn)) - 0.5 * sqrt(1. +
					    4.*(l*l*sigmaIn*sigmaIn)*deltaIn/ \
                        (AIn*AIn*sin(thetaIn)*sin(thetaIn))) - 
                        2.*A_SPIN*R_MIN*l/AIn);
		    } else {
			    lnOfh = 1.;
            }

            if (lnOfh <0. || r < R_MIN) {

                prim[j][i][RHO] = 1e-7*RHO_MIN;
                prim[j][i][UU] = 1e-7*U_MIN;
                prim[j][i][U1] = 0.;
                prim[j][i][U2] = 0.;
                prim[j][i][U3] = 0.;

            } else {

                prim[j][i][RHO] = pow((exp(lnOfh) - 1.)*
                    (ADIABATIC_INDEX-1.)/(KAPPA*ADIABATIC_INDEX),
                    1./(ADIABATIC_INDEX-1.));

                if (prim[j][i][RHO] > rhoMax)
                    rhoMax = prim[j][i][RHO];

                PetscRandomGetValue(randNumGen, &randNum);
                prim[j][i][UU] = KAPPA*pow(prim[j][i][RHO], ADIABATIC_INDEX)/\
                              (ADIABATIC_INDEX-1.)*(1. + 4e-2*(randNum-0.5));

                if (prim[j][i][UU] > uMax && r > R_MIN)
                    uMax = prim[j][i][UU];

                uPhiTmp = sqrt((-1. + sqrt(1. + 4.*l*l*(sigma*sigma*delta/\
                               (A*A*sin(theta)*sin(theta))))) / 2.);
			
                uconBL[1] = 0.;
                uconBL[2] = 0.;
                uconBL[3] = 2.*A_SPIN*r*sqrt(1. + uPhiTmp*uPhiTmp)/ \
                            sqrt(A*sigma*delta) + 
                            sqrt(sigma/A)*uPhiTmp/sin(theta);
                
                // Transform uconBoyerLinquist to uconKerrSchild and set to prim

                for (int ii=0; ii<NDIM; ii++)
                    for (int jj=0; jj<NDIM; jj++) {
                        gcovBL[ii][jj] = 0.;
                        gconBL[ii][jj] = 0.;
                        transform[ii][jj] = 0.;
                    }

                DD = 1. - 2./r + A_SPIN*A_SPIN/(r*r);
                mu = 1 + A_SPIN*A_SPIN*cos(theta)*cos(theta)/(r*r);

                gcovBL[0][0] = -(1. - 2./(r*mu));
                gcovBL[0][3] = -2.*A_SPIN*sin(theta)*sin(theta)/(r*mu);
                gcovBL[3][0] = gcovBL[0][3];
                gcovBL[1][1] = mu/DD;
                gcovBL[2][2] = r*r*mu;
                gcovBL[3][3] = r*r*sin(theta)*sin(theta)*\
                               (1. + A_SPIN*A_SPIN/(r*r) + 
                                2.*A_SPIN*A_SPIN*sin(theta)*sin(theta)/\
                                (r*r*r*mu));

                gconBL[0][0] = -1. -2.*(1 + A_SPIN*A_SPIN/(r*r))/(r*DD*mu);
                gconBL[0][3] = -2.*A_SPIN/(r*r*r*DD*mu);
                gconBL[3][0] = gconBL[0][3];
                gconBL[1][1] = DD/mu;
                gconBL[2][2] = 1./(r*r*mu);
                gconBL[3][3] = (1. - 2./(r*mu))/(r*r*sin(theta)*sin(theta)*DD);

                transform[0][0] = 1.;
                transform[1][1] = 1.;
                transform[2][2] = 1.;
                transform[3][3] = 1.;
                transform[0][1] = 2.*r/(r*r* - 2.*r + A_SPIN*A_SPIN);
                transform[3][1] = A_SPIN/(r*r - 2.*r + A_SPIN*A_SPIN);

                AA = gcovBL[0][0];
                BB = 2.*(gcovBL[0][1]*uconBL[1] + 
                         gcovBL[0][2]*uconBL[2] +
                         gcovBL[0][3]*uconBL[3]);
                CC = 1. + gcovBL[1][1]*uconBL[1]*uconBL[1] +
                          gcovBL[2][2]*uconBL[2]*uconBL[2] +
                          gcovBL[3][3]*uconBL[3]*uconBL[3] +
                      2.*(gcovBL[1][2]*uconBL[1]*uconBL[2] +
                          gcovBL[1][3]*uconBL[1]*uconBL[3] +
                          gcovBL[2][3]*uconBL[2]*uconBL[3]);

                discriminent = BB*BB - 4.*AA*CC;
                uconBL[0] = -(BB + sqrt(discriminent))/(2.*AA);

                uconKS[0] = transform[0][0]*uconBL[0] + 
                            transform[0][1]*uconBL[1] +
                            transform[0][2]*uconBL[2] +
                            transform[0][3]*uconBL[3];

                uconKS[1] = transform[1][0]*uconBL[0] + 
                            transform[1][1]*uconBL[1] +
                            transform[1][2]*uconBL[2] +
                            transform[1][3]*uconBL[3];

                uconKS[2] = transform[2][0]*uconBL[0] + 
                            transform[2][1]*uconBL[1] +
                            transform[2][2]*uconBL[2] +
                            transform[2][3]*uconBL[3];

                uconKS[3] = transform[3][0]*uconBL[0] + 
                            transform[3][1]*uconBL[1] +
                            transform[3][2]*uconBL[2] +
                            transform[3][3]*uconBL[3];

                PetscScalar rFactor, hFactor;
                rFactor = r - R0;
                hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*X2);
                uconKSPrime[0] = uconKS[0];
                uconKSPrime[1] = uconKS[1]*(1./rFactor);
                uconKSPrime[2] = uconKS[2]*(1./hFactor);
                uconKSPrime[3] = uconKS[3];

                PetscScalar gcov[NDIM][NDIM], gcon[NDIM][NDIM];
                PetscScalar gdet, alpha;
                gCovCalc(gcov, X1, X2);
                gDetCalc(&gdet, gcov);
                gConCalc(gcon, gcov, gdet);
                alphaCalc(&alpha, gcon);

                prim[j][i][U1] = uconKSPrime[1] + 
                                 alpha*alpha*uconKSPrime[0]*gcon[0][1];
                prim[j][i][U2] = uconKSPrime[2] + 
                                 alpha*alpha*uconKSPrime[0]*gcon[0][2];
                prim[j][i][U3] = uconKSPrime[3] +
                                 alpha*alpha*uconKSPrime[0]*gcon[0][3];

            }

            prim[j][i][B1] = 0.;
            prim[j][i][B2] = 0.;
            prim[j][i][B3] = 0.;

        }


    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        
            prim[j][i][RHO] = prim[j][i][RHO]/rhoMax;
            prim[j][i][UU] = prim[j][i][UU]/rhoMax;

            AVector[j][i] = 0.;
        }
    }

    uMax = uMax/rhoMax;

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            REAL gcov[NDIM][NDIM];
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloorInit(vars, gamma, X1, X2);

            for (int var=0; var<DOF; var++)
                prim[j][i][var] = vars[var];
        }
    }

    for (int j=X2Start-1; j<X2Start+X2Size+1; j++) {
        for (int i=X1Start-1; i<X1Start+X1Size+1; i++) {        
            rhoAv = 0.25*(prim[j][i][RHO] + prim[j][i-1][RHO] +
                          prim[j-1][i][RHO] + prim[j-1][i-1][RHO]);
            
            q = rhoAv - 0.2;

            if (q > 0.)
                AVector[j][i] = q;

        }
    }

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
            
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            PetscScalar gcov[NDIM][NDIM], gcon[NDIM][NDIM];
            PetscScalar g, gdet, alpha;
            gCovCalc(gcov, X1, X2);
            gDetCalc(&gdet, gcov);
            gConCalc(gcon, gcov, gdet);
            alphaCalc(&alpha, gcon);
            g = sqrt(-gdet);

            prim[j][i][B1] = -(AVector[j][i] - AVector[j+1][i] +
                               AVector[j][i+1] - AVector[j+1][i+1])/\
                              (2.*DX2*g);

            prim[j][i][B2] = (AVector[j][i] + AVector[j+1][i] -
                            AVector[j][i+1] - AVector[j+1][i+1])/\
                            (2.*DX1*g);

            prim[j][i][B3] = 0.;

            PetscScalar gamma, ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
            PetscScalar bsqr, var[DOF];

            for (int n=0; n<DOF; n++)
                var[n] = prim[j][i][n];

            gammaCalc(&gamma, var, gcov);
            uconCalc(ucon, gamma, alpha, var, gcon);
            covFromCon(ucov, ucon, gcov);
            bconCalc(bcon, var, ucon, ucov);
            covFromCon(bcov, bcon, gcov);
            bSqrCalc(&bsqr, bcon, bcov);

            if (bsqr > bsqrMax)
                bsqrMax = bsqr;

        }
    }
    betaActual = (ADIABATIC_INDEX-1.)*uMax/(0.5*bsqrMax);
    norm = sqrt(betaActual/BETA);

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
            prim[j][i][B1] = prim[j][i][B1]*norm;
            prim[j][i][B2] = prim[j][i][B2]*norm;
        }
    }

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            REAL gcov[NDIM][NDIM];
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloorInit(vars, gamma, X1, X2);

            prim[j][i][RHO] = (vars[RHO]);
            prim[j][i][UU] = (vars[UU]);

            for (int var=2; var<DOF; var++)
                prim[j][i][var] = vars[var];

        }
    }

    DMLocalToGlobalBegin(dmda, localPrim, INSERT_VALUES, Prim);
    DMLocalToGlobalEnd(dmda, localPrim, INSERT_VALUES, Prim);

    DMDAVecRestoreArrayDOF(dmda, localPrim, &prim);
    DMRestoreLocalVector(dmda, &localPrim);

    PetscRandomDestroy(&randNumGen);
}


