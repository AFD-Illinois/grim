#include "grim.h"

PetscErrorCode IFunction(TS ts, 
                         PetscScalar t,
                         Vec Xvec, Vec dX_dt,
                         Vec F, void *ptr)
{
    AppCtx *ctx = (AppCtx*)ptr;

    DM da;

    TSGetDM(ts, &da);

    DMDAGetCorners(da, 
                   &ctx->xstart, &ctx->ystart, &ctx->zstart,
                   &ctx->xsize, &ctx->ysize, &ctx->zsize);
    

    switch (ctx->dim) {
        case 1:
        {
            PetscScalar **x, **f;
            PetscScalar **dx_dt;

            DMDAVecGetArrayDOF(da, Xvec, &x);
            DMDAVecGetArrayDOF(da, F, &f);
            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);

            ComputedU_dt1D(f, x, dx_dt, ctx);

            DMDAVecRestoreArrayDOF(da, Xvec, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);
            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);

            break;
        }
        case 2:
        {
            PetscScalar ***x, ***f;
            PetscScalar ***dx_dt;

            DMDAVecGetArrayDOF(da, Xvec, &x);
            DMDAVecGetArrayDOF(da, F, &f);
            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);

            ComputedU_dt2D(f, x, dx_dt, ctx);

            DMDAVecRestoreArrayDOF(da, Xvec, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);
            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);

            break;
        }
        case 3:
        {
            PetscScalar ****x, ****f;
            PetscScalar ****dx_dt;

            DMDAVecGetArrayDOF(da, Xvec, &x);
            DMDAVecGetArrayDOF(da, F, &f);
            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);

            ComputedU_dt3D(f, x, dx_dt, ctx);

            DMDAVecRestoreArrayDOF(da, Xvec, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);
            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);

            break;
        }
        default: SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim in IFunction\n");
    }

    return(0);
}
