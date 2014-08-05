struct consToPrimInverter
{
    SNES snes;
    DM dmdaSnes;
    Vec residualPetscVec;
    struct geometry *geom;
    struct fluidElem *elem;
    REAL conservedVars[DOF];
}

void consToPrimInverterInit
(
    struct consToPrimInverter *inverter
)
{
    SNESCreate(PETSC_COMM_WORLD, &inverter->snes);
    VecCreate(PETSC_COMM_WORLD, &inverter->residualPetscVec);
    VecSetSize(inverter->residualPetscVec, DOF, DOF);
}

void setConsToPrimInverterData(struct geometry *geom,
                               REAL conservedVars[DOF],
                               struct consToPrimInverter *inverter)
{
    inverter->geom = geom;
    inverter->conservedVars = conservedVars;
}

PetscErrorCode consToPrimInverterResidual
(
    SNES snes,
    Vec primVarsGuessPetscVec,
    Vec residualPetscVec,
    struct consToPrimInverter *inverter
)
{
    REAL *primVarsGues, *residual;
    REAL conservedVarsGuess[DOF];

    VecGetArrayRead(primVarsGuessPetscVec, &primVarsGues);
    VecGetArray(residualPetscVec, &residual);

    setFluidElement(primVarsGuess, inverter->geom, inverter->elem);

    computeFluxes(inverter->elem, inverter->geom, 0, conservedVarsGuess);

    for (int var=0; var<DOF; var++)
    {
        residual[var] = conservedVarsGuess[var] - inverter->conservedVars[var];
    }

}
