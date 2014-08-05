void RiemannSolver(REAL flux[DOF],
                   const REAL fluxL[DOF],
                   const REAL fluxR[DOF],
                   const REAL uL[DOF],
                   const REAL uR[DOF],
                   const REAL cminL,
                   const REAL cmaxL,
                   const REAL cminR,
                   const REAL cmaxR)
{
	REAL cmax = fabs(fmax(fmax(0., cmaxL), cmaxR));
	REAL cmin = fabs(fmax(fmax(0., -cminL), -cminR));
	REAL ctop = fmax(cmax, cmin);
    
    for (int var=0; var<DOF; var++) {
        flux[var] = 0.5*(fluxL[var] + fluxR[var] -
                         ctop*(uR[var] - uL[var]));
    }
}

