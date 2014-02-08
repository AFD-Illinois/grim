void RiemannSolver(const REAL* restrict fluxL,
                   const REAL* restrict fluxR,
                   const REAL* restrict uL,
                   const REAL* restrict uR,
                   REAL* restrict flux)
{
    for (int var=0; var<DOF; var++) {
        flux[var] = 0.5*(fluxL[var] + fluxR[var] +
                         (uL[var] - uR[var]));
    }
}
