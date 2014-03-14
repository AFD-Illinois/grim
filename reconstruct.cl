#define LEFT -1
#define RIGHT 1
#define DOWN -1
#define UP 1

REAL SlopeLim(REAL y1, REAL y2, REAL y3)
{
    REAL Dqm = 2. * (y2 - y1);
	REAL Dqp = 2. * (y3 - y2);
	REAL Dqc = 0.5 * (y3 - y1);
	REAL s = Dqm * Dqp;
	if (s <= 0.) {
		return 0.;
    }
	else {
		if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
			return (Dqm);
		else if (fabs(Dqp) < fabs(Dqc))
			return (Dqp);
		else
			return (Dqc);
	}


//	REAL Dqm = (y2 - y1);
//	REAL Dqp = (y3 - y2);
//	REAL s = Dqm * Dqp;
//	if (s <= 0.)
//	    return 0.;
//	else if (fabs(Dqm) < fabs(Dqp))
//		return Dqm;
//	else
//		return Dqp;
        

//	REAL Dqm = (y2 - y1);
//	REAL Dqp = (y3 - y2);
//	REAL s = Dqm * Dqp;
//	if (s <= 0.)
//	    return 0.;
//	else
//	    return (2. * s / (Dqm + Dqp));
}

void ReconstructX1(const __local REAL* restrict primTile, 
                   const int iTile, const int jTile,
                   REAL* restrict primEdge,
                   const int dir)
{
    REAL left, mid, right, slope;

//    left = exp(primTile[INDEX_LOCAL(iTile-1, jTile, RHO)]);
//    mid = exp(primTile[INDEX_LOCAL(iTile, jTile, RHO)]);
//    right = exp(primTile[INDEX_LOCAL(iTile+1, jTile, RHO)]);
//    
//    slope = SlopeLim(left, mid, right);
//
//    primEdge[RHO] = mid + dir*0.5*slope;
//
//
//    left = exp(primTile[INDEX_LOCAL(iTile-1, jTile, UU)]);
//    mid = exp(primTile[INDEX_LOCAL(iTile, jTile, UU)]);
//    right = exp(primTile[INDEX_LOCAL(iTile+1, jTile, UU)]);
//    
//    slope = SlopeLim(left, mid, right);
//
//    primEdge[UU] = mid + dir*0.5*slope;

    for (int var=0; var<DOF; var++) {
        left = primTile[INDEX_LOCAL(iTile-1, jTile, var)];
        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
        right = primTile[INDEX_LOCAL(iTile+1, jTile, var)];
        
        slope = SlopeLim(left, mid, right);

        primEdge[var] = mid + dir*0.5*slope;
    }
}

void ReconstructX2(const __local REAL* restrict primTile, 
                   const int iTile, const int jTile,
                   REAL* restrict primEdge,
                   const int dir)
{
    REAL left, mid, right, slope;

//    left = exp(primTile[INDEX_LOCAL(iTile, jTile-1, RHO)]);
//    mid = exp(primTile[INDEX_LOCAL(iTile, jTile, RHO)]);
//    right = exp(primTile[INDEX_LOCAL(iTile, jTile+1, RHO)]);
//    
//    slope = SlopeLim(left, mid, right);
//
//    primEdge[RHO] = mid + dir*0.5*slope;
//
//
//    left = exp(primTile[INDEX_LOCAL(iTile, jTile-1, UU)]);
//    mid = exp(primTile[INDEX_LOCAL(iTile, jTile, UU)]);
//    right = exp(primTile[INDEX_LOCAL(iTile, jTile+1, UU)]);
//    
//    slope = SlopeLim(left, mid, right);
//
//    primEdge[UU] = mid + dir*0.5*slope;

    for (int var=0; var<DOF; var++) {
        left = primTile[INDEX_LOCAL(iTile, jTile-1, var)];
        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
        right = primTile[INDEX_LOCAL(iTile, jTile+1, var)];
        
        slope = SlopeLim(left, mid, right);

        primEdge[var] = mid + dir*0.5*slope;
    }

}
