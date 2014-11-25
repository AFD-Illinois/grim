#include "reconstruct.h"

/* Taken from HARM */
void slopeLim(const REAL left[ARRAY_ARGS DOF],
              const REAL mid[ARRAY_ARGS DOF],
              const REAL right[ARRAY_ARGS DOF],
              REAL ans[ARRAY_ARGS DOF])
{
  for (int var=0; var<DOF; var++)
  {
    /* Monotonized Central slope limiter */
    REAL Dqm = 2. * (mid[var] - left[var]);
	  REAL Dqp = 2. * (right[var] - mid[var]);
	  REAL Dqc = 0.5 * (right[var] - left[var]);
	  REAL s = Dqm * Dqp;
	  if (s <= 0.) 
    {
		  ans[var] = 0.;
    }
	  else
    {
		  if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      {
			  ans[var] = Dqm;
      }
		  else if (fabs(Dqp) < fabs(Dqc))
      {
			  ans[var] = Dqp;
      }
		  else
      {
			  ans[var] = Dqc;
      }
	  }

  }

}

/* The following code taken from PLUTO. Needs to be cleaned up. */
#define MINMOD(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
/*! Return the maximum between two numbers. */
#define MAX(a,b)  ( (a) >= (b) ? (a) : (b) ) 

/*! Return the minimum between two numbers. */
#define MIN(a,b)  ( (a) <= (b) ? (a) : (b) ) 
/* ************************************************************* */
double Median (double a, double b, double c)
/*
 *
 *
 *
 *************************************************************** */
{
  return (a + MINMOD(b-a,c-a));
}
REAL MP5_Reconstruct(const REAL Fjm2,
                     const REAL Fjm1,
                     const REAL Fj,
                     const REAL Fjp1,
                     const REAL Fjp2,
                     const REAL dx)
{
  double f, d2, d2p, d2m; 
  double dMMm, dMMp;
  double scrh1,scrh2, fmin, fmax; 
  double fAV, fMD, fLC, fUL, fMP;
  static double alpha = 4.0, epsm = 1.e-12;

  f  = 2.0*Fjm2 - 13.0*Fjm1 + 47.0*Fj + 27.0*Fjp1 - 3.0*Fjp2;
  f /= 60.0;   

  fMP = Fj + MINMOD(Fjp1-Fj,alpha*(Fj-Fjm1));

  if ((f - Fj)*(f - fMP) <= epsm) return f;

  d2m = Fjm2 + Fj - 2.0*Fjm1;    /* -- Eq. (2.19) -- */
  d2  = Fjm1 + Fjp1 - 2.0*Fj;
  d2p = Fj + Fjp2 - 2.0*Fjp1;    /* -- Eq. (2.19) -- */

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  fUL = Fj + alpha*(Fj - Fjm1);   /* -- Eq. (2.8) -- */
  fAV = 0.5*(Fj + Fjp1);        /* -- Eq. (2.16) -- */
  fMD = fAV - 0.5*dMMp; /* -- Eq. (2.28) -- */
  fLC = 0.5*(3.0*Fj - Fjm1) + 4.0/3.0*dMMm;  /* -- Eq. (2.29) -- */

  scrh1 = MIN(Fj, Fjp1); scrh1 = MIN(scrh1, fMD);
  scrh2 = MIN(Fj, fUL);    scrh2 = MIN(scrh2, fLC);
  fmin  = MAX(scrh1, scrh2);  /* -- Eq. (2.24a) -- */

  scrh1 = MAX(Fj, Fjp1); scrh1 = MAX(scrh1, fMD);
  scrh2 = MAX(Fj, fUL);    scrh2 = MAX(scrh2, fLC);
  fmax  = MIN(scrh1, scrh2);  /* -- Eq. (2.24b) -- */

  f = Median(f, fmin, fmax); /* -- Eq. (2.26) -- */
  return f;
}
#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))
#define MY_SIGN(fval) ( ((fval) <0.) ? -1. : 1. )
void para(const REAL x1,
          const REAL x2,
          const REAL x3,
          const REAL x4,
          const REAL x5,
          REAL lout[ARRAY_ARGS 1], REAL rout[ARRAY_ARGS 1])
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

         y[0]=x1;
         y[1]=x2;
         y[2]=x3;
         y[3]=x4;
         y[4]=x5;

         /*CW1.7 */
         for(i=1 ; i<4 ; i++) {
               Dqm = 2. *(y[i]-y[i-1]);
               Dqp = 2. *(y[i+1]-y[i]);
               Dqc = 0.5 *(y[i+1]-y[i-1]);
               aDqm = fabs(Dqm) ;
               aDqp = fabs(Dqp) ;
               aDqc = fabs(Dqc) ;
               s = Dqm*Dqp;

               if (s <=0.) dq[i]=0.;       //CW1.8
               else dq[i]=MY_MIN(aDqc,MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
         }

         /* CW1.6 */
         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;

         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));

         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

         if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
         else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;

         *lout=l;   //a_L,j
         *rout=r;
}
/*
 * left and right state reconsitruction using WENO-5,
 * all numbers from WHAM paper
 *
 * using zone-centered value of n=5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 */

/* author: Monika Moscibrodzka */
/* modified and sign error corrected by CFG 12.28.13 */

void weno(const REAL x1,
          const REAL x2,
          const REAL x3,
          const REAL x4,
          const REAL x5,
          REAL rout[ARRAY_ARGS 1])
{

    /* based on shu scholarpedia article, eqs 1,2,3 */
    /* 3rd order interpolations */
    double vr[3];
    vr[0] =  (3./8.)*x1-(5./4.)*x2+(15./8.)*x3;
    vr[1] = (-1./8.)*x2+(3./4.)*x3+(3./8.)*x4;
    vr[2] =  (3./8.)*x3+(3./4.)*x4-(1./8.)*x5;

    /* smoothness indicators, from Tchekh et al. A18, equiv to Shu eq. 8 */
    double beta[3];
    /*
    beta[0] = (4*x1*x1)/3.-(19*x1*x2)/3.+(25*x2*x2)/3.+(11*x1*x3)/3.-(31*x2*x3)/3.+(10*x3*x3)/3.;
    beta[1]= (4*x2*x2)/3.-(13*x2*x3)/3.+(13*x3*x3)/3.+(5*x2*x4)/3.-(13*x3*x4)/3.+(4*x4*x4)/3.;
    beta[2] = (4*x5*x5)/3.-(19*x5*x4)/3.+(25*x4*x4)/3.+(11*x5*x3)/3.-(31*x4*x3)/3.+(10*x3*x3)/3.;
    */
    beta[0]=(13./12.)*pow(x1-2.*x2+x3,2)+(1./4.)*pow(x1-4.*x2+3.*x3,2);
    beta[1]=(13./12.)*pow(x2-2.*x3+x4,2)+(1./4.)*pow(x4-x2,2);
    beta[2]=(13./12.)*pow(x3-2.*x4+x5,2)+(1./4.)*pow(x5-4.*x4+3.*x3,2);
    
    /* nonlinear weights, after shu eq. 9 */
    double den,wtr[3],Wr,wr[3],eps;
    eps=1.e-26;

    den = eps+beta[0]; den *= den; wtr[0] = (1./16.)/den;
    den = eps+beta[1]; den *= den; wtr[1] = (5./8.)/den;
    den = eps+beta[2]; den *= den; wtr[2] = (5./16.)/den;
    Wr = wtr[0]+wtr[1]+wtr[2];
    wr[0] = wtr[0]/Wr ;
    wr[1] = wtr[1]/Wr ;
    wr[2] = wtr[2]/Wr ;

    *rout = vr[0]*wr[0]+vr[1]*wr[1]+vr[2]*wr[2];
}

void reconstruct(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                 const int dir,
                 const int iTile, const int jTile,
                 const int X1Start, const int X2Start,
                 const int X1Size, const int X2Size,
                 REAL primVarsLeft[ARRAY_ARGS TILE_SIZE],
                 REAL primVarsRight[ARRAY_ARGS TILE_SIZE])
{

  if (dir==X1)
  {
    LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);
    
      REAL slope[DOF];
  
      slopeLim(&primTile[INDEX_TILE_MINUS_ONE_X1(&zone, 0)],
               &primTile[INDEX_TILE(&zone, 0)],
               &primTile[INDEX_TILE_PLUS_ONE_X1(&zone, 0)],
               slope);

      for (int var=0; var<DOF; var++)
      {
        /* Left Edge */
        primVarsLeft[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE(&zone, var)] - 0.5*slope[var];

        /* Right Edge */
        primVarsRight[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE(&zone, var)] + 0.5*slope[var];
      }
    
//      for (int var=0; var<DOF; var++)
//      {
//        primVarsRight[INDEX_TILE(&zone, var)] = 
//          MP5_Reconstruct(
//              primTile[INDEX_TILE_MINUS_TWO_X1(&zone, var)],
//              primTile[INDEX_TILE_MINUS_ONE_X1(&zone, var)],
//              primTile[INDEX_TILE(&zone, var)],
//              primTile[INDEX_TILE_PLUS_ONE_X1(&zone, var)],
//              primTile[INDEX_TILE_PLUS_TWO_X1(&zone, var)], zone.dX1);
//
//        primVarsLeft[INDEX_TILE(&zone, var)] = 
//          MP5_Reconstruct(
//              primTile[INDEX_TILE_PLUS_TWO_X1(&zone, var)],
//              primTile[INDEX_TILE_PLUS_ONE_X1(&zone, var)],
//              primTile[INDEX_TILE(&zone, var)],
//              primTile[INDEX_TILE_MINUS_ONE_X1(&zone, var)],
//              primTile[INDEX_TILE_MINUS_TWO_X1(&zone, var)], zone.dX1);
//
//      }

    }

  }
  else if (dir==X2)
  {
    LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);
    
      REAL slope[DOF];
  
      slopeLim(&primTile[INDEX_TILE_MINUS_ONE_X2(&zone, 0)],
               &primTile[INDEX_TILE(&zone, 0)],
               &primTile[INDEX_TILE_PLUS_ONE_X2(&zone, 0)],
               slope);

      for (int var=0; var<DOF; var++)
      {
        /* Left Edge */
        primVarsLeft[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE(&zone, var)] - 0.5*slope[var];

        /* Right Edge */
        primVarsRight[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE(&zone, var)] + 0.5*slope[var];
      }

//      for (int var=0; var<DOF; var++)
//      {
//        primVarsRight[INDEX_TILE(&zone, var)] = 
//          MP5_Reconstruct(
//              primTile[INDEX_TILE_MINUS_TWO_X2(&zone, var)],
//              primTile[INDEX_TILE_MINUS_ONE_X2(&zone, var)],
//              primTile[INDEX_TILE(&zone, var)],
//              primTile[INDEX_TILE_PLUS_ONE_X2(&zone, var)],
//              primTile[INDEX_TILE_PLUS_TWO_X2(&zone, var)], zone.dX2);
//
//        primVarsLeft[INDEX_TILE(&zone, var)] = 
//          MP5_Reconstruct(
//              primTile[INDEX_TILE_PLUS_TWO_X2(&zone, var)],
//              primTile[INDEX_TILE_PLUS_ONE_X2(&zone, var)],
//              primTile[INDEX_TILE(&zone, var)],
//              primTile[INDEX_TILE_MINUS_ONE_X2(&zone, var)],
//              primTile[INDEX_TILE_MINUS_TWO_X2(&zone, var)], zone.dX2);
//
//      }

    }

  }

}
