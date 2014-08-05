#include "grim.h"

PetscErrorCode SetCoordinates(DM dmda)
{
    DM coordDM;
    Vec coordVec;

    int X1Start, X2Start;
    int X1Size, X2Size;

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    DMDASetUniformCoordinates(dmda,
                              0., 1., 0., 1., 0., 1.);
    DMGetCoordinateDM(dmda, &coordDM);
    DMGetCoordinatesLocal(dmda, &coordVec);

    PetscScalar X1, X2, r, theta;
    DMDACoor2d **coord;
    
    DMDAVecGetArray(coordDM, coordVec, &coord);

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
    
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            BLCoords(&r, &theta, X1, X2);

            coord[j][i].x = r*PetscSinScalar(theta);
            coord[j][i].y = r*PetscCosScalar(theta);
        }
    }
    DMDAVecRestoreArray(coordDM, coordVec, &coord);

    DMSetCoordinatesLocal(dmda, coordVec);

    return (0);
}

