#include "grim.h"

PetscErrorCode RHSFunction(TS ts,
                           PetscScalar t, 
                           Vec Xvec, Vec F,
                           void *ptr)
{
    AppCtx *ctx = (AppCtx*)ptr;
    
    DM da;

    Vec localX;

    TSGetDM(ts, &da);

    DMGetLocalVector(da, &localX);

    DMDAGetCorners(da, 
                   &ctx->xstart, &ctx->ystart, &ctx->zstart,
                   &ctx->xsize, &ctx->ysize, &ctx->zsize);

    DMGlobalToLocalBegin(da, Xvec, INSERT_VALUES, localX);
    DMGlobalToLocalEnd(da, Xvec, INSERT_VALUES, localX);


    switch (ctx->dim) {
        case 1:
        {
            PetscInt i, var;
            PetscScalar **x, **f;
            PetscScalar **prim_lx, **prim_rx;
            PetscScalar **flux_lx, **flux_rx;
            PetscScalar **flux_x;
            PetscScalar **U_lx, **U_rx;

            DMDAVecGetArrayDOF(da, localX, &x);
            DMDAVecGetArrayDOF(da, F, &f);

            DMDAVecGetArrayDOF(da, localX, &x);
            DMDAVecGetArrayDOF(da, F, &f);

            DMDAVecGetArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecGetArrayDOF(da, ctx->prim_r[X], &prim_rx);

            DMDAVecGetArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecGetArrayDOF(da, ctx->flux_r[X], &flux_rx);

            DMDAVecGetArrayDOF(da, ctx->flux[X], &flux_x);

            DMDAVecGetArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecGetArrayDOF(da, ctx->U_r[X], &U_rx);

            Reconstruct1D(prim_lx, prim_rx, x, ctx);

            ComputeFlux1D(flux_lx, prim_lx, 1, ctx);
            ComputeFlux1D(flux_rx, prim_rx, 1, ctx);

            ComputeFlux1D(U_lx, prim_lx, 0, ctx);
            ComputeFlux1D(U_rx, prim_rx, 0, ctx);

            RiemannSolver1D(flux_x, 
                            flux_lx, flux_rx,
                            U_lx, U_rx, ctx);

            for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++)
                for (var=0; var<dof; var++)
                    f[i][var] = -(flux_x[i+1][var] - flux_x[i][var])/ctx->dx;

            DMDAVecRestoreArrayDOF(da, localX, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);

            DMDAVecRestoreArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecRestoreArrayDOF(da, ctx->prim_r[X], &prim_rx);

            DMDAVecRestoreArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecRestoreArrayDOF(da, ctx->flux_r[X], &flux_rx);

            DMDAVecRestoreArrayDOF(da, ctx->flux[X], &flux_x);

            DMDAVecRestoreArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecRestoreArrayDOF(da, ctx->U_r[X], &U_rx);
            break;
        }
        case 2:
        {
            PetscInt i, j, var;
            PetscScalar ***x, ***f;
            PetscScalar ***prim_lx, ***prim_rx;
            PetscScalar ***prim_ly, ***prim_ry;
            PetscScalar ***flux_lx, ***flux_rx;
            PetscScalar ***flux_ly, ***flux_ry;
            PetscScalar ***flux_x, ***flux_y;
            PetscScalar ***U_lx, ***U_rx;
            PetscScalar ***U_ly, ***U_ry;

            DMDAVecGetArrayDOF(da, localX, &x);
            DMDAVecGetArrayDOF(da, F, &f);

            DMDAVecGetArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecGetArrayDOF(da, ctx->prim_r[X], &prim_rx);
            DMDAVecGetArrayDOF(da, ctx->prim_l[Y], &prim_ly);
            DMDAVecGetArrayDOF(da, ctx->prim_r[Y], &prim_ry);

            DMDAVecGetArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecGetArrayDOF(da, ctx->flux_r[X], &flux_rx);
            DMDAVecGetArrayDOF(da, ctx->flux_l[Y], &flux_ly);
            DMDAVecGetArrayDOF(da, ctx->flux_r[Y], &flux_ry);

            DMDAVecGetArrayDOF(da, ctx->flux[X], &flux_x);
            DMDAVecGetArrayDOF(da, ctx->flux[Y], &flux_y);

            DMDAVecGetArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecGetArrayDOF(da, ctx->U_r[X], &U_rx);
            DMDAVecGetArrayDOF(da, ctx->U_l[Y], &U_ly);
            DMDAVecGetArrayDOF(da, ctx->U_r[Y], &U_ry);

            Reconstruct2D(prim_lx, prim_rx,
                          prim_ly, prim_ry,
                          x, ctx);

            ComputeFlux2D(flux_lx, prim_lx, 1, ctx);
            ComputeFlux2D(flux_rx, prim_rx, 1, ctx);
            ComputeFlux2D(flux_ly, prim_ly, 2, ctx);
            ComputeFlux2D(flux_ry, prim_ry, 2, ctx);

            ComputeFlux2D(U_lx, prim_lx, 0, ctx);
            ComputeFlux2D(U_rx, prim_rx, 0, ctx);
            ComputeFlux2D(U_ly, prim_ly, 0, ctx);
            ComputeFlux2D(U_ry, prim_ry, 0, ctx);

            RiemannSolver2D(flux_x, flux_y,
                            flux_lx, flux_rx,
                            flux_ly, flux_ry,
                            U_lx, U_rx, 
                            U_ly, U_ry, ctx);

            for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
                for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++)
                    for (var=0; var<dof; var++) {
                        f[j][i][var] = -(flux_x[j][i+1][var] - 
                                        flux_x[j][i][var])/ctx->dx +
                                       -(flux_y[j+1][i][var] -
                                        flux_y[j][i][var])/ctx->dy;

                    }


            DMDAVecRestoreArrayDOF(da, localX, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);

            DMDAVecRestoreArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecRestoreArrayDOF(da, ctx->prim_r[X], &prim_rx);
            DMDAVecRestoreArrayDOF(da, ctx->prim_l[Y], &prim_ly);
            DMDAVecRestoreArrayDOF(da, ctx->prim_r[Y], &prim_ry);

            DMDAVecRestoreArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecRestoreArrayDOF(da, ctx->flux_r[X], &flux_rx);
            DMDAVecRestoreArrayDOF(da, ctx->flux_l[Y], &flux_ly);
            DMDAVecRestoreArrayDOF(da, ctx->flux_r[Y], &flux_ry);

            DMDAVecRestoreArrayDOF(da, ctx->flux[X], &flux_x);
            DMDAVecRestoreArrayDOF(da, ctx->flux[Y], &flux_y);

            DMDAVecRestoreArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecRestoreArrayDOF(da, ctx->U_r[X], &U_rx);
            DMDAVecRestoreArrayDOF(da, ctx->U_l[Y], &U_ly);
            DMDAVecRestoreArrayDOF(da, ctx->U_r[Y], &U_ry);

            break;
        }
        case 3:
        {
            PetscInt i, j, k, var;
            PetscScalar ****x, ****f;
            PetscScalar ****prim_lx, ****prim_rx;
            PetscScalar ****prim_ly, ****prim_ry;
            PetscScalar ****prim_lz, ****prim_rz;
            PetscScalar ****flux_lx, ****flux_rx;
            PetscScalar ****flux_ly, ****flux_ry;
            PetscScalar ****flux_lz, ****flux_rz;
            PetscScalar ****U_lx, ****U_rx;
            PetscScalar ****U_ly, ****U_ry;
            PetscScalar ****U_lz, ****U_rz;
            PetscScalar ****flux_x, ****flux_y, ****flux_z;

            DMDAVecGetArrayDOF(da, localX, &x);
            DMDAVecGetArrayDOF(da, F, &f);

            DMDAVecGetArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecGetArrayDOF(da, ctx->prim_r[X], &prim_rx);
            DMDAVecGetArrayDOF(da, ctx->prim_l[Y], &prim_ly);
            DMDAVecGetArrayDOF(da, ctx->prim_r[Y], &prim_ry);
            DMDAVecGetArrayDOF(da, ctx->prim_l[Z], &prim_lz);
            DMDAVecGetArrayDOF(da, ctx->prim_r[Z], &prim_rz);

            DMDAVecGetArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecGetArrayDOF(da, ctx->flux_r[X], &flux_rx);
            DMDAVecGetArrayDOF(da, ctx->flux_l[Y], &flux_ly);
            DMDAVecGetArrayDOF(da, ctx->flux_r[Y], &flux_ry);
            DMDAVecGetArrayDOF(da, ctx->flux_l[Z], &flux_lz);
            DMDAVecGetArrayDOF(da, ctx->flux_r[Z], &flux_rz);

            DMDAVecGetArrayDOF(da, ctx->flux[X], &flux_x);
            DMDAVecGetArrayDOF(da, ctx->flux[Y], &flux_y);
            DMDAVecGetArrayDOF(da, ctx->flux[Z], &flux_z);

            DMDAVecGetArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecGetArrayDOF(da, ctx->U_r[X], &U_rx);
            DMDAVecGetArrayDOF(da, ctx->U_l[Y], &U_ly);
            DMDAVecGetArrayDOF(da, ctx->U_r[Y], &U_ry);
            DMDAVecGetArrayDOF(da, ctx->U_l[Z], &U_lz);
            DMDAVecGetArrayDOF(da, ctx->U_r[Z], &U_rz);

            Reconstruct3D(prim_lx, prim_rx,
                          prim_ly, prim_ry,
                          prim_lz, prim_rz,
                          x, ctx);

            ComputeFlux3D(flux_lx, prim_lx, 1, ctx);
            ComputeFlux3D(flux_rx, prim_rx, 1, ctx);
            ComputeFlux3D(flux_ly, prim_ly, 2, ctx);
            ComputeFlux3D(flux_ry, prim_ry, 2, ctx);
            ComputeFlux3D(flux_lz, prim_lz, 3, ctx);
            ComputeFlux3D(flux_rz, prim_rz, 3, ctx);

            ComputeFlux3D(U_lx, prim_lx, 0, ctx);
            ComputeFlux3D(U_rx, prim_rx, 0, ctx);
            ComputeFlux3D(U_ly, prim_ly, 0, ctx);
            ComputeFlux3D(U_ry, prim_ry, 0, ctx);
            ComputeFlux3D(U_lz, prim_lz, 0, ctx);
            ComputeFlux3D(U_rz, prim_rz, 0, ctx);

            RiemannSolver3D(flux_x, flux_y, flux_z,
                            flux_lx, flux_rx,
                            flux_ly, flux_ry,
                            flux_lz, flux_rz,
                            U_lx, U_rx, 
                            U_ly, U_ry,
                            U_lz, U_rz, ctx);

            for (k=ctx->zstart; k<ctx->zstart + ctx->zsize; k++)
                for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
                    for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++)
                        for (var=0; var<dof; var++) {
                            f[k][j][i][var] = -(flux_x[k][j][i+1][var] - 
                                            flux_x[k][j][i][var])/ctx->dx +
                                           
                                           -(flux_y[k][j+1][i][var] -
                                            flux_y[k][j][i][var])/ctx->dy +

                                           -(flux_z[k+1][j][i][var] -
                                            flux_z[k][j][i][var])/ctx->dz;

                        }

            DMDAVecRestoreArrayDOF(da, localX, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);

            DMDAVecRestoreArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecRestoreArrayDOF(da, ctx->prim_r[X], &prim_rx);
            DMDAVecRestoreArrayDOF(da, ctx->prim_l[Y], &prim_ly);
            DMDAVecRestoreArrayDOF(da, ctx->prim_r[Y], &prim_ry);
            DMDAVecRestoreArrayDOF(da, ctx->prim_l[Z], &prim_lz);
            DMDAVecRestoreArrayDOF(da, ctx->prim_r[Z], &prim_rz);

            DMDAVecRestoreArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecRestoreArrayDOF(da, ctx->flux_r[X], &flux_rx);
            DMDAVecRestoreArrayDOF(da, ctx->flux_l[Y], &flux_ly);
            DMDAVecRestoreArrayDOF(da, ctx->flux_r[Y], &flux_ry);
            DMDAVecRestoreArrayDOF(da, ctx->flux_l[Z], &flux_lz);
            DMDAVecRestoreArrayDOF(da, ctx->flux_r[Z], &flux_rz);

            DMDAVecRestoreArrayDOF(da, ctx->flux[X], &flux_x);
            DMDAVecRestoreArrayDOF(da, ctx->flux[Y], &flux_y);
            DMDAVecRestoreArrayDOF(da, ctx->flux[Z], &flux_z);

            DMDAVecRestoreArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecRestoreArrayDOF(da, ctx->U_r[X], &U_rx);
            DMDAVecRestoreArrayDOF(da, ctx->U_l[Y], &U_ly);
            DMDAVecRestoreArrayDOF(da, ctx->U_r[Y], &U_ry);
            DMDAVecRestoreArrayDOF(da, ctx->U_l[Z], &U_lz);
            DMDAVecRestoreArrayDOF(da, ctx->U_r[Z], &U_rz);

            break;
        }

        default: SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim in RHSFunction\n");

    }


    DMRestoreLocalVector(da, &localX);
    return(0);
}
