//void ComputedU_dt(const REAL* restrict var,
//                  const REAL* restrict dvar_dt,
//                  REAL* restrict dU_dt,
//                  REAL X1, REAL X2)
//{
//    REAL g, gDet, gamma, dgamma_dt, alpha, tmp1, tmp2, dtmp1_dt, dtmp2_dt;
//
//    gDet = GDET(X1,X2); g = sqrt(-gDet);
//    gamma = GAMMA; dgamma_dt = DGAMMA_DT; alpha = ALPHA;
//    tmp1 = P + var[RHO] + var[UU] + BSQR;
//    tmp2 = P + 0.5*BSQR;
//    dtmp1_dt = DP_DT + dvar_dt[RHO] + dvar_dt[UU] + DBSQR_DT;
//    dtmp2_dt = DP_DT + 0.5*DBSQR_DT;
//
//
//    dU_dt[RHO] = g*(dvar_dt[RHO]*UCON0 + var[RHO]*DUCON0_DT);
//
//    dU_dt[UU] = g*(dtmp1_dt*UCON0*UCOV0 +
//                   tmp1*(DUCON0_DT*UCOV0 + UCON0*DUCOV0_DT) +
//                   dtmp2_dt - DBCON0_DT*BCOV0 - BCON0*DBCOV0_DT);
//
//    dU_dt[U1] = g*(dtmp1_dt*UCON0*UCOV1 +
//                   tmp1*(DUCON0_DT*UCOV1 + UCON0*DUCOV1_DT) -
//                   DBCON0_DT*BCOV1 - BCON0*DBCOV1_DT);
//
//    dU_dt[U2] = g*(dtmp1_dt*UCON0*UCOV2 +
//                   tmp1*(DUCON0_DT*UCOV2 + UCON0*DUCOV2_DT) -
//                   DBCON0_DT*BCOV2 - BCON0*DBCOV2_DT);
//
//    dU_dt[U3] = g*(dtmp1_dt*UCON0*UCOV3 +
//                   tmp1*(DUCON0_DT*UCOV3 + UCON0*DUCOV3_DT) -
//                   DBCON0_DT*BCOV3 - BCON0*DBCOV3_DT);
//
//    dU_dt[B1] = g*dvar_dt[B1];
//    dU_dt[B2] = g*dvar_dt[B2];
//    dU_dt[B3] = g*dvar_dt[B3];
//
//}
//
//void addSources(const REAL* restrict var,
//                REAL* restrict dU_dt,
//                REAL X1, REAL X2)
//{
//    REAL g, gDet, gamma, alpha, tmp1, tmp2;
//    
//    gDet = GDET(X1,X2); g = sqrt(-gDet);
//    gamma = GAMMA; alpha = ALPHA;
//    tmp1 = P + var[RHO] + var[UU] + BSQR;
//    tmp2 = P + 0.5*BSQR;
//
//    REAL gh[NDIM][NDIM], gl[NDIM][NDIM], gcon[NDIM][NDIM];
//    REAL conntmp[NDIM][NDIM][NDIM], conn[NDIM][NDIM][NDIM];
//    REAL Xh[NDIM], Xl[NDIM], X[NDIM], MHD[NDIM][NDIM];
//    int i, j, k, l;
//
//    X[0] = 0.; X[1] = X1;  X[2] = X2; X[3] = 0.;
//
//	for (k = 0; k < NDIM; k++) {
//		for (l = 0; l < NDIM; l++)
//			Xh[l] = X[l];
//		for (l = 0; l < NDIM; l++)
//			Xl[l] = X[l];
//		Xh[k] += EPS;
//		Xl[k] -= EPS;
//        
//        gh[0][0] = GCOV00(Xh[1], Xh[2]);
//        gh[0][1] = GCOV01(Xh[1], Xh[2]);
//        gh[0][2] = GCOV02(Xh[1], Xh[2]);
//        gh[0][3] = GCOV03(Xh[1], Xh[2]);
//        gh[1][0] = GCOV10(Xh[1], Xh[2]);
//        gh[1][1] = GCOV11(Xh[1], Xh[2]);
//        gh[1][2] = GCOV12(Xh[1], Xh[2]);
//        gh[1][3] = GCOV13(Xh[1], Xh[2]);
//        gh[2][0] = GCOV20(Xh[1], Xh[2]);
//        gh[2][1] = GCOV21(Xh[1], Xh[2]);
//        gh[2][2] = GCOV22(Xh[1], Xh[2]);
//        gh[2][3] = GCOV23(Xh[1], Xh[2]);
//        gh[3][0] = GCOV30(Xh[1], Xh[2]);
//        gh[3][1] = GCOV31(Xh[1], Xh[2]);
//        gh[3][2] = GCOV32(Xh[1], Xh[2]);
//        gh[3][3] = GCOV33(Xh[1], Xh[2]);
//
//        gl[0][0] = GCOV00(Xl[1], Xl[2]);
//        gl[0][1] = GCOV01(Xl[1], Xl[2]);
//        gl[0][2] = GCOV02(Xl[1], Xl[2]);
//        gl[0][3] = GCOV03(Xl[1], Xl[2]);
//        gl[1][0] = GCOV10(Xl[1], Xl[2]);
//        gl[1][1] = GCOV11(Xl[1], Xl[2]);
//        gl[1][2] = GCOV12(Xl[1], Xl[2]);
//        gl[1][3] = GCOV13(Xl[1], Xl[2]);
//        gl[2][0] = GCOV20(Xl[1], Xl[2]);
//        gl[2][1] = GCOV21(Xl[1], Xl[2]);
//        gl[2][2] = GCOV22(Xl[1], Xl[2]);
//        gl[2][3] = GCOV23(Xl[1], Xl[2]);
//        gl[3][0] = GCOV30(Xl[1], Xl[2]);
//        gl[3][1] = GCOV31(Xl[1], Xl[2]);
//        gl[3][2] = GCOV32(Xl[1], Xl[2]);
//        gl[3][3] = GCOV33(Xl[1], Xl[2]);
//
//
//        gcon[0][0] = GCON00(X[1], X[2], gDet);
//        gcon[0][1] = GCON01(X[1], X[2], gDet);
//        gcon[0][2] = GCON02(X[1], X[2], gDet);
//        gcon[0][3] = GCON03(X[1], X[2], gDet);
//        gcon[1][0] = GCON10(X[1], X[2], gDet);
//        gcon[1][1] = GCON11(X[1], X[2], gDet);
//        gcon[1][2] = GCON12(X[1], X[2], gDet);
//        gcon[1][3] = GCON13(X[1], X[2], gDet);
//        gcon[2][0] = GCON20(X[1], X[2], gDet);
//        gcon[2][1] = GCON21(X[1], X[2], gDet);
//        gcon[2][2] = GCON22(X[1], X[2], gDet);
//        gcon[2][3] = GCON23(X[1], X[2], gDet);
//        gcon[3][0] = GCON30(X[1], X[2], gDet);
//        gcon[3][1] = GCON31(X[1], X[2], gDet);
//        gcon[3][2] = GCON32(X[1], X[2], gDet);
//        gcon[3][3] = GCON33(X[1], X[2], gDet);
//
//
//		for (i = 0; i < NDIM; i++)
//			for (j = 0; j < NDIM; j++)
//				conn[i][j][k] =
//				    (gh[i][j] - gl[i][j]) / (Xh[k] -
//							     Xl[k]);
//	}
//
//	/* now rearrange to find \Gamma_{ijk} */
//	for (i = 0; i < NDIM; i++)
//		for (j = 0; j < NDIM; j++)
//			for (k = 0; k < NDIM; k++)
//				conntmp[i][j][k] =
//				    0.5 * (conn[j][i][k] + conn[k][i][j] -
//					   conn[k][j][i]);
//
//	/* finally, raise index */
//	for (i = 0; i < NDIM; i++)
//		for (j = 0; j < NDIM; j++)
//			for (k = 0; k < NDIM; k++) {
//				conn[i][j][k] = 0.;
//				for (l = 0; l < NDIM; l++)
//					conn[i][j][k] +=
//					    gcon[i][l] *
//					    conntmp[l][j][k];
//			}
//    
//    MHD[0][0] = MHD0_0;
//    MHD[0][1] = MHD0_1;
//    MHD[0][2] = MHD0_2;
//    MHD[0][3] = MHD0_3;
//    MHD[1][0] = MHD1_0;
//    MHD[1][1] = MHD1_1;
//    MHD[1][2] = MHD1_2;
//    MHD[1][3] = MHD1_3;
//    MHD[2][0] = MHD2_0;
//    MHD[2][1] = MHD2_1;
//    MHD[2][2] = MHD2_2;
//    MHD[2][3] = MHD2_3;
//    MHD[3][0] = MHD3_0;
//    MHD[3][1] = MHD3_1;
//    MHD[3][2] = MHD3_2;
//    MHD[3][3] = MHD3_3;
//
//    for (j=0; j<NDIM; j++)
//        for (k=0; k<NDIM; k++) {
//            dU_dt[UU] = dU_dt[UU] - g*(MHD[j][k]*conn[k][0][j]);
//            dU_dt[U1] = dU_dt[U1] - g*(MHD[j][k]*conn[k][1][j]);
//            dU_dt[U2] = dU_dt[U2] - g*(MHD[j][k]*conn[k][2][j]);
//            dU_dt[U3] = dU_dt[U3] - g*(MHD[j][k]*conn[k][3][j]);
//    
//        }
//
//
//}
//
//void ComputeFluxAndUX1(const REAL* restrict var,
//                       REAL* restrict flux,
//                       REAL* restrict U,
//                       REAL X1, REAL X2)
//{
//    REAL g, gDet, gamma, alpha, tmp1, tmp2;
//    
//    gDet = GDET(X1,X2); g = sqrt(-gDet);
//    gamma = GAMMA; alpha = ALPHA;
//    tmp1 = P + var[RHO] + var[UU] + BSQR;
//    tmp2 = P + 0.5*BSQR;
//
//    flux[RHO]  = g*var[RHO]*UCON1;
//
//    flux[UU] = g*MHD1_0;
//    flux[U1] = g*MHD1_1;
//    flux[U2] = g*MHD1_2;
//    flux[U3] = g*MHD1_3;
//
//    flux[B1] = 0;
//    flux[B2] = g*(BCON2*UCON1 - BCON1*UCON2);
//    flux[B3] = g*(BCON3*UCON1 - BCON1*UCON3);
//
//    U[RHO]  = g*var[RHO]*UCON0;
//
//    U[UU] = g*MHD0_0;
//    U[U1] = g*MHD0_1;
//    U[U2] = g*MHD0_2;
//    U[U3] = g*MHD0_3;
//
//    U[B1] = g*(BCON1*UCON0 - BCON0*UCON1);
//    U[B2] = g*(BCON2*UCON0 - BCON0*UCON2);
//    U[B3] = g*(BCON3*UCON0 - BCON0*UCON3);
//
//}
//
//void ComputeFluxAndUX2(const REAL* restrict var,
//                       REAL* restrict flux,
//                       REAL* restrict U,
//                       REAL X1, REAL X2)
//{
//    REAL g, gDet, gamma, alpha, tmp1, tmp2;
//
//    gDet = GDET(X1,X2); g = sqrt(-gDet);
//    gamma = GAMMA; alpha = ALPHA;
//    tmp1 = P + var[RHO] + var[UU] + BSQR;
//    tmp2 = P + 0.5*BSQR;
//
//    flux[RHO]  = g*var[RHO]*UCON2;
//
//    flux[UU] = g*MHD2_0;
//    flux[U1] = g*MHD2_1;
//    flux[U2] = g*MHD2_2;
//    flux[U3] = g*MHD2_3;
//
//    flux[B1] = g*(BCON1*UCON2 - BCON2*UCON1);
//    flux[B2] = 0;
//    flux[B3] = g*(BCON3*UCON2 - BCON2*UCON3);
//
//    U[RHO]  = g*var[RHO]*UCON0;
//
//    U[UU] = g*MHD0_0;
//    U[U1] = g*MHD0_1;
//    U[U2] = g*MHD0_2;
//    U[U3] = g*MHD0_3;
//
//    U[B1] = g*(BCON1*UCON0 - BCON0*UCON1);
//    U[B2] = g*(BCON2*UCON0 - BCON0*UCON2);
//    U[B3] = g*(BCON3*UCON0 - BCON0*UCON3);
//
//}
//
