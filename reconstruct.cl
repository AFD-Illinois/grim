#define LEFT -1
#define RIGHT 1

REAL SlopeLim(REAL y1, REAL y2, REAL y3)
{
    REAL tmp1 = y2 - y1;
    REAL tmp2 = y3 - y2;
    REAL s = 2.*tmp1*tmp2;

    if (s<=0)
        return 0.;
    else
        return s/(tmp1 + tmp2);
}

void ReconstructX1(const __local REAL* restrict primTile, 
                   const int iTile, const int jTile,
                   REAL* restrict primEdge,
                   const int dir)
{
    REAL left, mid, right, slope;
    for (int var=0; var<DOF; var++) {
        left = primTile[INDEX_LOCAL(iTile-1, jTile, var)];
        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
        right = primTile[INDEX_LOCAL(iTile+1, jTile, var)];
        
        slope = SlopeLim(left, mid, right);

        primEdge[var] = mid + dir*slope;
    }
}

void ReconstructX2(const __local REAL* restrict primTile, 
                   const int iTile, const int jTile,
                   REAL* restrict primEdge,
                   const int dir)
{
    REAL left, mid, right, slope;
    for (int var=0; var<DOF; var++) {
        left = primTile[INDEX_LOCAL(iTile, jTile-1, var)];
        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
        right = primTile[INDEX_LOCAL(iTile, jTile+1, var)];
        
        slope = SlopeLim(left, mid, right);

        primEdge[var] = mid + dir*slope;
    }

}
