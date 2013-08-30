#include "grim.h"
#include <petscviewerhdf5.h>

enum {rho, u, ux, uy, uz, bx, by, bz, dof};
enum {X, Y, Z}; //Directions;

typedef enum {BDRY_PERIODIC, BDRY_OUTFLOW} bdry_condn;
static const char *const bdry_types[] =
{"PERIODIC", "OUTFLOW", "bdry_condn", "BDRY_", 0};

typedef enum {LIM_MC, LIM_VANL, LIM_MINM} limiter;
static const char *const limiter_types[] =
{"MC", "VANL", "MINM", "LIM_", 0};

typedef struct {
    PetscInt dim, Nx, Ny, Nz;
    PetscScalar dx, dy, dz;
    PassiveScalar xmin, xmax, ymin, ymax, zmin, zmax;
    PetscInt nghost;
    PassiveScalar dt, cfl;
    bdry_condn bdry_x, bdry_y, bdry_z;
    limiter lim;

    PetscInt xstart, xsize, ystart, ysize, zstart, zsize;
    PetscInt xstart_work_array, xsize_work_array;
    PetscInt ystart_work_array, ysize_work_array;
    PetscInt zstart_work_array, zsize_work_array;

    Vec *prim_l, *prim_r;
    Vec *flux_l, *flux_r;
    Vec *U_l, *U_r;
    Vec *flux, U;
} AppCtx;


PETSC_EXTERN PetscErrorCode InitialCondition(DM da, Vec X, void *ptr);
PETSC_EXTERN PetscErrorCode ResFunction(TS ts, PetscScalar t, 
                                        Vec Xvec, Vec dX_dt, Vec F, void *ptr);

PETSC_EXTERN PetscScalar slope_lim(PetscScalar y1, PetscScalar y2, 
                                   PetscScalar y3, AppCtx *ctx);

PETSC_EXTERN PetscErrorCode Reconstruct1D(PetscScalar **prim_lx, 
                                       PetscScalar **prim_rx,
                                       PetscScalar **x, AppCtx *ctx);

PETSC_EXTERN PetscErrorCode Reconstruct2D(PetscScalar ***prim_lx, 
                                       PetscScalar ***prim_rx,
                                       PetscScalar ***prim_ly,
                                       PetscScalar ***prim_ry,
                                       PetscScalar ***x, AppCtx *ctx);

PETSC_EXTERN PetscErrorCode Reconstruct3D(PetscScalar ****prim_lx, 
                                       PetscScalar ****prim_rx,
                                       PetscScalar ****prim_ly,
                                       PetscScalar ****prim_ry,
                                       PetscScalar ****prim_lz,
                                       PetscScalar ****prim_rz,
                                       PetscScalar ****x, AppCtx *ctx);

PETSC_EXTERN PetscErrorCode ComputeFlux1D(PetscScalar **flux, 
                                          PetscScalar **x, PetscInt dir, 
                                          AppCtx *ctx);

PETSC_EXTERN PetscErrorCode ComputeFlux2D(PetscScalar ***flux, 
                                          PetscScalar ***x, PetscInt dir, 
                                          AppCtx *ctx);

PETSC_EXTERN PetscErrorCode ComputeFlux3D(PetscScalar ****flux, 
                                          PetscScalar ****x, PetscInt dir, 
                                          AppCtx *ctx);


PETSC_EXTERN PetscErrorCode ComputedU_dt1D(PetscScalar **dU_dt, 
                                           PetscScalar **x, 
                                           PetscScalar **dx_dt,
                                           AppCtx *ctx);

PETSC_EXTERN PetscErrorCode ComputedU_dt2D(PetscScalar ***dU_dt, 
                                           PetscScalar ***x, 
                                           PetscScalar ***dx_dt,
                                           AppCtx *ctx);

PETSC_EXTERN PetscErrorCode ComputedU_dt3D(PetscScalar ****dU_dt, 
                                           PetscScalar ****x, 
                                           PetscScalar ****dx_dt,
                                           AppCtx *ctx);

PETSC_EXTERN PetscErrorCode RiemannSolver1D(PetscScalar **flux_x, 
                                            PetscScalar **flux_lx, 
                                            PetscScalar **flux_rx,
                                            PetscScalar **U_lx,
                                            PetscScalar **U_rx, AppCtx *ctx);

PETSC_EXTERN PetscErrorCode RiemannSolver2D(PetscScalar ***flux_x, 
                                            PetscScalar ***flux_y,
                                            PetscScalar ***flux_lx, 
                                            PetscScalar ***flux_rx,
                                            PetscScalar ***flux_ly, 
                                            PetscScalar ***flux_ry,
                                            PetscScalar ***U_lx, 
                                            PetscScalar ***U_rx,
                                            PetscScalar ***U_ly, 
                                            PetscScalar ***U_ry, AppCtx *ctx);

PETSC_EXTERN PetscErrorCode RiemannSolver3D(PetscScalar ****flux_x, 
                                            PetscScalar ****flux_y,
                                            PetscScalar ****flux_z,
                                            PetscScalar ****flux_lx, 
                                            PetscScalar ****flux_rx,
                                            PetscScalar ****flux_ly, 
                                            PetscScalar ****flux_ry,
                                            PetscScalar ****flux_lz, 
                                            PetscScalar ****flux_rz,
                                            PetscScalar ****U_lx, 
                                            PetscScalar ****U_rx,
                                            PetscScalar ****U_ly, 
                                            PetscScalar ****U_ry,
                                            PetscScalar ****U_lz, 
                                            PetscScalar ****U_rz, AppCtx *ctx);

#include "IFunction.cpp"
#include "RHSFunction.cpp"
int main(int argc, char **argv)
{
    AppCtx ctx; 
    TS ts;
    SNES snes;
    Vec soln;
    DM da;
    PetscErrorCode ierr;

    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    // Defaults
    ctx.dim = 1;
    ctx.xmin = 0.; ctx.xmax = 1.;
    ctx.ymin = 0.; ctx.ymax = 1.;
    ctx.zmin = 0.; ctx.zmax = 1.;
    ctx.cfl = 0.5;
    ctx.nghost = 2;
    ctx.bdry_x = BDRY_PERIODIC;
    ctx.bdry_y = BDRY_PERIODIC;
    ctx.bdry_z = BDRY_PERIODIC;
    ctx.lim = LIM_MC;
    // Initialization
    ctx.xstart = 0; ctx.xsize = 0;
    ctx.ystart = 0; ctx.ysize = 0;
    ctx.zstart = 0; ctx.zsize = 0;
    ctx.xstart_work_array= 0; ctx.xsize_work_array= 0;
    ctx.ystart_work_array= 0; ctx.ysize_work_array= 0;
    ctx.zstart_work_array= 0; ctx.zsize_work_array= 0;
    ctx.prim_l = PETSC_NULL; ctx.prim_r = PETSC_NULL;
    ctx.flux_l = PETSC_NULL; ctx.flux_r = PETSC_NULL;
    ctx.U_l = PETSC_NULL; ctx.U_r = PETSC_NULL;
    ctx.flux = PETSC_NULL; ctx.U = PETSC_NULL;

    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "GRIM options","");
    {
        PetscOptionsInt("-dim", "dim", "", ctx.dim, &ctx.dim, NULL);
        PetscOptionsInt("-nghost", "Number of ghost zones", "", ctx.nghost,
                        &ctx.nghost, NULL);
        PetscOptionsScalar("-xmin", "x0", "", ctx.xmin, &ctx.xmin, NULL);
        PetscOptionsScalar("-xmax", "x1", "", ctx.xmax, &ctx.xmax, NULL);
        PetscOptionsScalar("-ymin", "y0", "", ctx.ymin, &ctx.ymin, NULL);
        PetscOptionsScalar("-ymax", "y1", "", ctx.ymax, &ctx.ymax, NULL);
        PetscOptionsScalar("-zmin", "z0", "", ctx.zmin, &ctx.zmin, NULL);
        PetscOptionsScalar("-zmax", "z1", "", ctx.zmax, &ctx.zmax, NULL);
        PetscOptionsEnum("-bdry_x", "boundary in x", "", bdry_types, (PetscEnum)ctx.bdry_x,
                           (PetscEnum*)&ctx.bdry_x, NULL);
        PetscOptionsEnum("-bdry_y", "boundary in y", "", bdry_types, (PetscEnum)ctx.bdry_y,
                           (PetscEnum*)&ctx.bdry_y, NULL);
        PetscOptionsEnum("-bdry_z", "boundary in z", "", bdry_types, (PetscEnum)ctx.bdry_z,
                           (PetscEnum*)&ctx.bdry_z, NULL);
        PetscOptionsEnum("-limiter", "Slope limiter", "", limiter_types,
                         (PetscEnum)ctx.lim, (PetscEnum*)&ctx.lim, NULL);
        PetscOptionsScalar("-dt", "initial dt", "", ctx.dt, &ctx.dt, NULL);
    }
    PetscOptionsEnd();

    switch (ctx.dim) {
        case 1: 
        {
            ierr = DMDACreate1d(PETSC_COMM_WORLD,
                                DMDA_BOUNDARY_PERIODIC, -30, 
                                dof, ctx.nghost, PETSC_NULL, &da); CHKERRQ(ierr);

            DMDAGetInfo(da, 
                        PETSC_IGNORE, &ctx.Nx, PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
            ctx.dx = (ctx.xmax - ctx.xmin)/(PetscScalar)ctx.Nx;

            DMDASetUniformCoordinates(da, ctx.xmin + ctx.dx/2., ctx.xmax - ctx.dx/2.,
                                      0., 0., 
                                      0., 0.);
        
         
            break;

        }

        case 2:
        {
            ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                             DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                             DMDA_STENCIL_STAR,
                             -30, -30,
                             PETSC_DECIDE, PETSC_DECIDE,
                             dof, ctx.nghost, PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);

            DMDAGetInfo(da, 
                        PETSC_IGNORE, &ctx.Nx, &ctx.Ny, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
            ctx.dx = (ctx.xmax - ctx.xmin)/(PetscScalar)ctx.Nx;
            ctx.dy = (ctx.ymax - ctx.ymin)/(PetscScalar)ctx.Ny;

            DMDASetUniformCoordinates(da, ctx.xmin + ctx.dx/2., ctx.xmax - ctx.dx/2.,
                                          ctx.ymin + ctx.dy/2., ctx.ymax - ctx.dy/2., 
                                          0., 0.);
            break;
        }

        case 3: 
        {
            ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                             DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                             DMDA_BOUNDARY_PERIODIC,
                             DMDA_STENCIL_STAR,
                             -30, -30, -30,
                             PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                             dof, ctx.nghost, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da);
                             CHKERRQ(ierr);
            DMDAGetInfo(da, 
                        PETSC_IGNORE, &ctx.Nx, &ctx.Ny, &ctx.Nz, 
                        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, 
                        PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
            ctx.dx = (ctx.xmax - ctx.xmin)/(PetscScalar)ctx.Nx;
            ctx.dy = (ctx.ymax - ctx.ymin)/(PetscScalar)ctx.Ny;
            ctx.dz = (ctx.zmax - ctx.zmin)/(PetscScalar)ctx.Nz;

            DMDASetUniformCoordinates(da, ctx.xmin + ctx.dx/2., ctx.xmax - ctx.dx/2.,
                                          ctx.ymin + ctx.dy/2., ctx.ymax - ctx.dy/2., 
                                          ctx.zmin + ctx.dz/2., ctx.zmax - ctx.dz/2.);
            break;
        }
        default: SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim in creating DMDA\n");
    }
    DMCreateGlobalVector(da, &soln);

    DMDAGetGhostCorners(da, 
                   &ctx.xstart_work_array, &ctx.ystart_work_array, 
                   &ctx.zstart_work_array, &ctx.xsize_work_array,
                   &ctx.ysize_work_array, &ctx.zsize_work_array);

    DMCreateLocalVector(da, &ctx.U);

    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.prim_l);
    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.prim_r);

    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.U_l);
    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.U_r);

    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.flux_l);
    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.flux_r);

    VecDuplicateVecs(ctx.U, ctx.dim, &ctx.flux);

    DMDASetFieldName(da, rho, "rho");
    DMDASetFieldName(da, u, "u");
    DMDASetFieldName(da, ux, "ux");
    DMDASetFieldName(da, uy, "uy");
    DMDASetFieldName(da, uz, "uz");
    DMDASetFieldName(da, bx, "bx");
    DMDASetFieldName(da, by, "by");
    DMDASetFieldName(da, bz, "bz");

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetDM(ts, da);
    TSSetIFunction(ts, PETSC_NULL, ResFunction, &ctx);
//    TSSetIFunction(ts, PETSC_NULL, IFunction, &ctx);
//    TSSetRHSFunction(ts, PETSC_NULL, RHSFunction, &ctx);

    InitialCondition(da, soln, &ctx);
    TSSetSolution(ts, soln);
    TSSetType(ts, TSTHETA);
    TSSetFromOptions(ts);

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD,"initialcondition.h5", FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) soln, "var");
    VecView(soln, viewer);
    PetscViewerDestroy(&viewer);
    
    TSSolve(ts, soln);

    PetscViewerHDF5Open(PETSC_COMM_WORLD,"soln.h5", FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) soln, "var");
    VecView(soln, viewer);
    PetscViewerDestroy(&viewer);



    VecDestroy(&soln);
    VecDestroy(&ctx.U);

    VecDestroyVecs(ctx.dim, &ctx.prim_l);
    VecDestroyVecs(ctx.dim, &ctx.prim_r);

    VecDestroyVecs(ctx.dim, &ctx.U_l);
    VecDestroyVecs(ctx.dim, &ctx.U_r);

    VecDestroyVecs(ctx.dim, &ctx.flux_l);
    VecDestroyVecs(ctx.dim, &ctx.flux_r);

    VecDestroyVecs(ctx.dim, &ctx.flux);

    TSDestroy(&ts);
    DMDestroy(&da);

    PetscFinalize();
    return(0);
}

PetscErrorCode InitialCondition(DM da, Vec Xvec, void *ptr)
{
    AppCtx *ctx = (AppCtx*)ptr;

    DMDAGetCorners(da, 
                   &ctx->xstart, &ctx->ystart, &ctx->zstart,
                   &ctx->xsize, &ctx->ysize, &ctx->zsize);
                   
    switch (ctx->dim) {
        case 1:
        {
            PetscInt i;
            PetscScalar xcoord, r;

            PetscScalar **x;
            DMDAVecGetArrayDOF(da, Xvec, &x);
            
            for (i=ctx->xstart; i<ctx->xstart+ctx->xsize; i++) {
                xcoord = ctx->xmin + ctx->dx/2. + i*ctx->dx;
                r = xcoord - 0.5;
                x[i][rho] = 1. + 0.1*PetscSinScalar(2.*PETSC_PI*xcoord);
                x[i][u] = 0.01;
                x[i][ux] = 0.;
                x[i][uy] = 0.;
                x[i][uz] = 0.;
                x[i][bx] = 0.;
                x[i][by] = 0.;
                x[i][bz] = 0.;
            }

            DMDAVecRestoreArrayDOF(da, Xvec, &x);

            break;
        }

        case 2:
        {
            PetscInt i, j;
            PetscScalar xcoord, ycoord, r;
            PetscScalar vx, vy, gamma;

            PetscScalar ***x;
            DMDAVecGetArrayDOF(da, Xvec, &x);

            for (j=ctx->ystart; j<ctx->ystart+ctx->ysize; j++) {
                ycoord = ctx->ymin + ctx->dy/2. + j*ctx->dy;
                for (i=ctx->xstart; i<ctx->xstart+ctx->xsize; i++) {
                    xcoord = ctx->xmin + ctx->dx/2. + i*ctx->dx;
                    r = PetscSqrtScalar((xcoord-.5)*(xcoord-.5) + 
                                        (ycoord-.5)*(ycoord-.5));


//                    x[j][i][rho] = 1. + PetscExpScalar(-r*r/(0.1*0.1));
                    x[j][i][rho] = 25./(36.*M_PI);
                    x[j][i][u] = 5./(12.*M_PI*(5./3 - 1.));
                    vx = -0.5*PetscSinScalar(2*M_PI*ycoord);
                    vy = 0.5*PetscSinScalar(2*M_PI*xcoord);
                    gamma = 1./PetscSqrtScalar(1 - vx*vx - vy*vy);
                    x[j][i][ux] = gamma*vx;
                    x[j][i][uy] = gamma*vy;
                    x[j][i][uz] = 0.0;
                    x[j][i][bx] =
                        -1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(2*M_PI*ycoord);
                    x[j][i][by] =
                        1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(4*M_PI*xcoord);
                    x[j][i][bz] = 0.0;
                }
            }

            DMDAVecRestoreArrayDOF(da, Xvec, &x);

            break;
        }
        case 3:
        {
            PetscInt i, j, k;
            PetscScalar xcoord, ycoord, zcoord, r;

            PetscScalar ****x;
            DMDAVecGetArrayDOF(da, Xvec, &x);
            
            for (k=ctx->zstart; k<ctx->zstart+ctx->zsize; k++) {
                zcoord = ctx->zmin + ctx->dz/2. + k*ctx->dz;

                for (j=ctx->ystart; j<ctx->ystart+ctx->ysize; j++) {
                    ycoord = ctx->ymin + ctx->dy/2. + j*ctx->dy;
                
                    for (i=ctx->xstart; i<ctx->xstart+ctx->xsize; i++) {
                        xcoord = ctx->xmin + ctx->dx/2. + i*ctx->dx;
                        r = PetscSqrtScalar((xcoord-.5)*(xcoord-.5) + 
                                            (ycoord-.5)*(ycoord-.5) + 
                                            (zcoord-.5)*(zcoord-.5));

                        x[k][j][i][rho] = 1. + 0.1*PetscExpScalar(-r*r/0.01);
                        x[k][j][i][u] = 0.01;
                        x[k][j][i][ux] = 0.0;
                        x[k][j][i][uy] = 0.0;
                        x[k][j][i][uz] = 0.0;
                        x[k][j][i][bx] = 0.0;
                        x[k][j][i][by] = 0.0;
                        x[k][j][i][bz] = 0.0;
                    }
                }
            }

            DMDAVecRestoreArrayDOF(da, Xvec, &x);

            break;

        }

        default: SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim in initial conditions\n");

    }


    return(0);
}

PETSC_STATIC_INLINE PetscScalar delta(PetscInt up, PetscInt down)
{
    if (up!=down)
        return 0.;
    else return 1;
}

PetscErrorCode ComputeFlux1D(PetscScalar **flux, PetscScalar **x, 
                             PetscInt dir, AppCtx *ctx)
{
    PetscInt i;
    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;

    for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++) {
        P = (4./3 - 1)*x[i][u];

        gamma = PetscSqrtScalar(1. + x[i][ux]*x[i][ux] + x[i][uy]*x[i][uy] +
                                     x[i][uz]*x[i][uz]);
        u_con[0] = gamma;
        u_con[1] = x[i][ux];
        u_con[2] = x[i][uy];
        u_con[3] = x[i][uz];

        u_cov[0] = -gamma;
        u_cov[1] = x[i][ux];
        u_cov[2] = x[i][uy];
        u_cov[3] = x[i][uz];

        b_con[0] = x[i][bx]*u_con[1] + x[i][by]*u_con[2] + x[i][bz]*u_con[3];
        b_con[1] = (x[i][bx] + b_con[0]*u_con[1])/u_con[0];
        b_con[2] = (x[i][by] + b_con[0]*u_con[2])/u_con[0];
        b_con[3] = (x[i][bz] + b_con[0]*u_con[3])/u_con[0];

        b_cov[0] = -b_con[0];
        b_cov[1] = b_con[1];
        b_cov[2] = b_con[2];
        b_cov[3] = b_con[3];

        b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
                 b_con[3]*b_con[3];

        flux[i][rho] = x[i][rho]*u_con[dir];
        
        flux[i][bx] = b_con[1]*u_con[dir] - b_con[dir]*u_con[1];
        flux[i][by] = b_con[2]*u_con[dir] - b_con[dir]*u_con[2];
        flux[i][bz] = b_con[3]*u_con[dir] - b_con[dir]*u_con[3];

        flux[i][u] = (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*u_cov[0] + 
                     (P + 0.5*b_sqr)*delta(dir, 0) - b_con[dir]*b_cov[0];

        flux[i][ux] = (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*u_cov[1] + 
                     (P + 0.5*b_sqr)*delta(dir, 1) - b_con[dir]*b_cov[1];

        flux[i][uy] = (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*u_cov[2] + 
                     (P + 0.5*b_sqr)*delta(dir, 2) - b_con[dir]*b_cov[2];

        flux[i][uz] = (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*u_cov[3] + 
                     (P + 0.5*b_sqr)*delta(dir, 3) - b_con[dir]*b_cov[3];
    }

    return(0);

}

PetscErrorCode ComputeFlux2D(PetscScalar ***flux, PetscScalar ***x,
                             PetscInt dir, AppCtx *ctx)
{
    PetscInt i, j;
    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;

    for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++) 
        for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++) {
            P = (5./3 - 1)*x[j][i][u];

            gamma = PetscSqrtScalar(1. + x[j][i][ux]*x[j][i][ux] + x[j][i][uy]*x[j][i][uy] +
                                         x[j][i][uz]*x[j][i][uz]);
            u_con[0] = gamma;
            u_con[1] = x[j][i][ux];
            u_con[2] = x[j][i][uy];
            u_con[3] = x[j][i][uz];

            u_cov[0] = -gamma;
            u_cov[1] = x[j][i][ux];
            u_cov[2] = x[j][i][uy];
            u_cov[3] = x[j][i][uz];

            b_con[0] = x[j][i][bx]*u_con[1] + x[j][i][by]*u_con[2] + x[j][i][bz]*u_con[3];
            b_con[1] = (x[j][i][bx] + b_con[0]*u_con[1])/u_con[0];
            b_con[2] = (x[j][i][by] + b_con[0]*u_con[2])/u_con[0];
            b_con[3] = (x[j][i][bz] + b_con[0]*u_con[3])/u_con[0];

            b_cov[0] = -b_con[0];
            b_cov[1] = b_con[1];
            b_cov[2] = b_con[2];
            b_cov[3] = b_con[3];

            b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
                     b_con[3]*b_con[3];

            flux[j][i][rho] = x[j][i][rho]*u_con[dir];
        
            flux[j][i][bx] = b_con[1]*u_con[dir] - b_con[dir]*u_con[1];
            flux[j][i][by] = b_con[2]*u_con[dir] - b_con[dir]*u_con[2];
            flux[j][i][bz] = b_con[3]*u_con[dir] - b_con[dir]*u_con[3];

            flux[j][i][u] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[0] + 
                            (P + 0.5*b_sqr)*delta(dir, 0) - b_con[dir]*b_cov[0];

            flux[j][i][ux] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[1] + 
                             (P + 0.5*b_sqr)*delta(dir, 1) - b_con[dir]*b_cov[1];

            flux[j][i][uy] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[2] + 
                             (P + 0.5*b_sqr)*delta(dir, 2) - b_con[dir]*b_cov[2];

            flux[j][i][uz] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[3] + 
                             (P + 0.5*b_sqr)*delta(dir, 3) - b_con[dir]*b_cov[3];
        }

    return(0);

}

PetscErrorCode ComputeFlux3D(PetscScalar ****flux, PetscScalar ****x,
                             PetscInt dir, AppCtx *ctx)
{
    PetscInt i, j, k;
    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;

    for (k=ctx->zstart-1; k<ctx->zstart + ctx->zsize+1; k++)
        for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++)
            for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++) {
                P = (4./3 - 1)*x[k][j][i][u];

                gamma = PetscSqrtScalar(1. + x[k][j][i][ux]*x[k][j][i][ux] + x[k][j][i][uy]*x[k][j][i][uy] +
                                             x[k][j][i][uz]*x[k][j][i][uz]);
                u_con[0] = gamma;
                u_con[1] = x[k][j][i][ux];
                u_con[2] = x[k][j][i][uy];
                u_con[3] = x[k][j][i][uz];

                u_cov[0] = -gamma;
                u_cov[1] = x[k][j][i][ux];
                u_cov[2] = x[k][j][i][uy];
                u_cov[3] = x[k][j][i][uz];

                b_con[0] = x[k][j][i][bx]*u_con[1] + x[k][j][i][by]*u_con[2] + x[k][j][i][bz]*u_con[3];
                b_con[1] = (x[k][j][i][bx] + b_con[0]*u_con[1])/u_con[0];
                b_con[2] = (x[k][j][i][by] + b_con[0]*u_con[2])/u_con[0];
                b_con[3] = (x[k][j][i][bz] + b_con[0]*u_con[3])/u_con[0];

                b_cov[0] = -b_con[0];
                b_cov[1] = b_con[1];
                b_cov[2] = b_con[2];
                b_cov[3] = b_con[3];

                b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
                         b_con[3]*b_con[3];

                flux[k][j][i][rho] = x[k][j][i][rho]*u_con[dir];
        
                flux[k][j][i][bx] = b_con[1]*u_con[dir] - b_con[dir]*u_con[1];
                flux[k][j][i][by] = b_con[2]*u_con[dir] - b_con[dir]*u_con[2];
                flux[k][j][i][bz] = b_con[3]*u_con[dir] - b_con[dir]*u_con[3];

                flux[k][j][i][u] = (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*u_cov[0] + 
                                (P + 0.5*b_sqr)*delta(dir, 0) - b_con[dir]*b_cov[0];

                flux[k][j][i][ux] = (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*u_cov[1] + 
                                 (P + 0.5*b_sqr)*delta(dir, 1) - b_con[dir]*b_cov[1];

                flux[k][j][i][uy] = (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*u_cov[2] + 
                                 (P + 0.5*b_sqr)*delta(dir, 2) - b_con[dir]*b_cov[2];

                flux[k][j][i][uz] = (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*u_cov[3] + 
                                 (P + 0.5*b_sqr)*delta(dir, 3) - b_con[dir]*b_cov[3];
            }

    return(0);

}

PetscScalar slope_lim(PetscScalar y1, PetscScalar y2, PetscScalar y3, 
                      AppCtx *ctx)
{
    PetscScalar Dqm, Dqp, Dqc, s;

    /* Monotonized Central or Woodward limiter */
	if (ctx->lim == LIM_MC) {
		Dqm = 2. * (y2 - y1);
		Dqp = 2. * (y3 - y2);
		Dqc = 0.5 * (y3 - y1);
		s = Dqm * Dqp;
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
	}
	/* van leer slope limiter */
	else if (ctx->lim == LIM_VANL) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else
			return (2. * s / (Dqm + Dqp));
	}

	/* minmod slope limiter (crude but robust) */
	else if (ctx->lim == LIM_MINM) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else if (fabs(Dqm) < fabs(Dqp))
			return Dqm;
		else
			return Dqp;
	}

    SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim slope limiter\n");
    
}

PetscErrorCode Reconstruct1D(PetscScalar **prim_lx, PetscScalar **prim_rx,
                             PetscScalar **x, AppCtx *ctx)
{
    if (ctx->dim!=1) 
        SETERRQ1(PETSC_COMM_SELF, 1, "1D Reconstruct called for %d dim\n", ctx->dim);

    PetscInt i, var;
    PetscScalar slope;

    for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++)
        for (var=0; var<dof; var++) {
            slope = slope_lim(x[i-1][var], x[i][var], x[i+1][var], ctx);
                    
            prim_lx[i][var] = x[i][var] - 0.5*slope;
            prim_rx[i][var] = x[i][var] + 0.5*slope;

        }

    return(0);
}

PetscErrorCode Reconstruct2D(PetscScalar ***prim_lx, PetscScalar ***prim_rx,
                             PetscScalar ***prim_ly, PetscScalar ***prim_ry,
                             PetscScalar ***x, AppCtx *ctx)
{
    if (ctx->dim!=2) 
        SETERRQ1(PETSC_COMM_SELF, 1, "2D Reconstruct called for %d dim\n", ctx->dim);

    PetscInt i, j, var;
    PetscScalar slope;

    for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++)
        for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++)
            for (var=0; var<dof; var++) {
                slope = slope_lim(x[j][i-1][var], x[j][i][var],
                                  x[j][i+1][var], ctx);

                prim_lx[j][i][var] = x[j][i][var] - 0.5*slope;
                prim_rx[j][i][var] = x[j][i][var] + 0.5*slope;

                slope = slope_lim(x[j-1][i][var], x[j][i][var],
                                  x[j+1][i][var], ctx);
                    
                prim_ly[j][i][var] = x[j][i][var] - 0.5*slope;
                prim_ry[j][i][var] = x[j][i][var] + 0.5*slope;

            }

    return(0);
}

PetscErrorCode Reconstruct3D(PetscScalar ****prim_lx, PetscScalar ****prim_rx,
                             PetscScalar ****prim_ly, PetscScalar ****prim_ry,
                             PetscScalar ****prim_lz, PetscScalar ****prim_rz,
                             PetscScalar ****x, AppCtx *ctx)
{
    if (ctx->dim!=3) 
        SETERRQ1(PETSC_COMM_SELF, 1, "3D Reconstruct called for %d dim\n", ctx->dim);

    PetscInt i, j, k, var;
    PetscScalar slope;

    for (k=ctx->zstart-1; k<ctx->zstart + ctx->zsize+1; k++)
        for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++)
            for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++)
                for (var=0; var<dof; var++) {
                    slope = slope_lim(x[k][j][i-1][var],
                                      x[k][j][i][var],
                                      x[k][j][i+1][var], ctx);
                    
                    prim_lx[k][j][i][var] = x[k][j][i][var] - 0.5*slope;
                    prim_rx[k][j][i][var] = x[k][j][i][var] + 0.5*slope;

                    slope = slope_lim(x[k][j-1][i][var],
                                      x[k][j][i][var], 
                                      x[k][j+1][i][var], ctx);
                    
                    prim_ly[k][j][i][var] = x[k][j][i][var] - 0.5*slope;
                    prim_ry[k][j][i][var] = x[k][j][i][var] + 0.5*slope;

                    slope = slope_lim(x[k-1][j][i][var],
                                      x[k][j][i][var], 
                                      x[k+1][j][i][var], ctx);
                    
                    prim_lz[k][j][i][var] = x[k][j][i][var] - 0.5*slope;
                    prim_rz[k][j][i][var] = x[k][j][i][var] + 0.5*slope;

                }
    return(0);

}

PetscErrorCode ComputedU_dt1D(PetscScalar **dU_dt, 
                              PetscScalar **x, PetscScalar **dx_dt,
                              AppCtx *ctx)
{
    PetscInt i, dir=0;
    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;
    PetscScalar dgamma_dt, du_cov_dt[4], du_con_dt[4];
    PetscScalar db_cov_dt[4], db_con_dt[4], db_sqr_dt, dP_dt;

    for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++) {
        P = (4./3 - 1)*x[i][u];
        dP_dt = (4./3 - 1)*dx_dt[i][u];

        gamma = PetscSqrtScalar(1. + x[i][ux]*x[i][ux] + x[i][uy]*x[i][uy] +
                                     x[i][uz]*x[i][uz]);
        dgamma_dt = (x[i][ux]*dx_dt[i][ux] + x[i][uy]*dx_dt[i][uy] +
                    x[i][uz]*dx_dt[i][uz])/gamma;

        u_con[0] = gamma;
        u_con[1] = x[i][ux];
        u_con[2] = x[i][uy];
        u_con[3] = x[i][uz];

        du_con_dt[0] = dgamma_dt;
        du_con_dt[1] = dx_dt[i][ux];
        du_con_dt[2] = dx_dt[i][uy];
        du_con_dt[3] = dx_dt[i][uz];

        u_cov[0] = -gamma;
        u_cov[1] = x[i][ux];
        u_cov[2] = x[i][uy];
        u_cov[3] = x[i][uz];

        du_cov_dt[0] = -dgamma_dt;
        du_cov_dt[1] = dx_dt[i][ux];
        du_cov_dt[2] = dx_dt[i][uy];
        du_cov_dt[3] = dx_dt[i][uz];

        b_con[0] = x[i][bx]*u_con[1] + x[i][by]*u_con[2] + x[i][bz]*u_con[3];
        b_con[1] = (x[i][bx] + b_con[0]*u_con[1])/u_con[0];
        b_con[2] = (x[i][by] + b_con[0]*u_con[2])/u_con[0];
        b_con[3] = (x[i][bz] + b_con[0]*u_con[3])/u_con[0];

        db_con_dt[0] = dx_dt[i][bx]*u_con[1] + dx_dt[i][by]*u_con[2] + 
                       dx_dt[i][bz]*u_con[3] + x[i][bx]*du_con_dt[1] + 
                       x[i][by]*du_con_dt[2] + x[i][bz]*du_con_dt[3];

        db_con_dt[1] = -(x[i][bx] + b_con[0]*u_con[1])*du_con_dt[0]/
                        (u_con[0]*u_con[0]) + 
                        (dx_dt[i][bx] + b_con[0]*du_con_dt[1] + 
                         db_con_dt[0]*u_con[1])/u_con[0];
        
        db_con_dt[2] = -(x[i][by] + b_con[0]*u_con[2])*du_con_dt[0]/
                        (u_con[0]*u_con[0]) + 
                        (dx_dt[i][by] + b_con[0]*du_con_dt[2] + 
                         db_con_dt[0]*u_con[2])/u_con[0];

        db_con_dt[3] = -(x[i][bz] + b_con[0]*u_con[3])*du_con_dt[0]/
                        (u_con[0]*u_con[0]) + 
                        (dx_dt[i][bz] + b_con[0]*du_con_dt[3] + 
                         db_con_dt[0]*u_con[3])/u_con[0];

        b_cov[0] = -b_con[0];
        b_cov[1] = b_con[1];
        b_cov[2] = b_con[2];
        b_cov[3] = b_con[3];

        db_cov_dt[0] = -db_con_dt[0];
        db_cov_dt[1] = db_con_dt[1];
        db_cov_dt[2] = db_con_dt[2];
        db_cov_dt[3] = db_con_dt[3];

        b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
                 b_con[3]*b_con[3];

        db_sqr_dt = 2.*(-b_con[0]*db_con_dt[0] + b_con[1]*db_con_dt[1] + 
                        b_con[2]*db_con_dt[2] + b_con[3]*db_con_dt[3]);

        dU_dt[i][rho] = dx_dt[i][rho]*u_con[dir] + x[i][rho]*du_con_dt[dir];

        dU_dt[i][bx] = dx_dt[i][bx];

        dU_dt[i][by] = dx_dt[i][by];

        dU_dt[i][bz] = dx_dt[i][bz];

        dU_dt[i][u] = (dP_dt + dx_dt[i][rho] + dx_dt[i][u] +
                       db_sqr_dt)*u_con[dir]*u_cov[0] + 
                      (P + x[i][rho] + x[i][u] + b_sqr)*du_con_dt[dir]*u_cov[0]
                      + (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*du_cov_dt[0]
                      + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 0)
                      - db_con_dt[dir]*b_cov[0] - b_con[dir]*db_cov_dt[0];

        dU_dt[i][ux] = (dP_dt + dx_dt[i][rho] + dx_dt[i][u] +
                       db_sqr_dt)*u_con[dir]*u_cov[1] + 
                      (P + x[i][rho] + x[i][u] + b_sqr)*du_con_dt[dir]*u_cov[1]
                      + (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*du_cov_dt[1]
                      + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 1)
                      - db_con_dt[dir]*b_cov[1] - b_con[dir]*db_cov_dt[1];
        
        dU_dt[i][uy] = (dP_dt + dx_dt[i][rho] + dx_dt[i][u] +
                       db_sqr_dt)*u_con[dir]*u_cov[2] + 
                      (P + x[i][rho] + x[i][u] + b_sqr)*du_con_dt[dir]*u_cov[2]
                      + (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*du_cov_dt[2]
                      + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 2)
                      - db_con_dt[dir]*b_cov[2] - b_con[dir]*db_cov_dt[2];

        dU_dt[i][uz] = (dP_dt + dx_dt[i][rho] + dx_dt[i][u] +
                       db_sqr_dt)*u_con[dir]*u_cov[3] + 
                      (P + x[i][rho] + x[i][u] + b_sqr)*du_con_dt[dir]*u_cov[3]
                      + (P + x[i][rho] + x[i][u] + b_sqr)*u_con[dir]*du_cov_dt[3]
                      + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 3)
                      - db_con_dt[dir]*b_cov[3] - b_con[dir]*db_cov_dt[3];

    }

    return(0);

}

PetscErrorCode ComputedU_dt2D(PetscScalar ***dU_dt, 
                              PetscScalar ***x, PetscScalar ***dx_dt,
                              AppCtx *ctx)
{
    PetscInt i, j, dir=0;
    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;
    PetscScalar dgamma_dt, du_cov_dt[4], du_con_dt[4];
    PetscScalar db_cov_dt[4], db_con_dt[4], db_sqr_dt, dP_dt;

    for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
        for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++) {
            P = (5./3 - 1)*x[j][i][u];
            dP_dt = (5./3 - 1)*dx_dt[j][i][u];

            gamma = PetscSqrtScalar(1. + x[j][i][ux]*x[j][i][ux] + x[j][i][uy]*x[j][i][uy] +
                                         x[j][i][uz]*x[j][i][uz]);
            dgamma_dt = (x[j][i][ux]*dx_dt[j][i][ux] + x[j][i][uy]*dx_dt[j][i][uy] +
                         x[j][i][uz]*dx_dt[j][i][uz])/gamma;

            u_con[0] = gamma;
            u_con[1] = x[j][i][ux];
            u_con[2] = x[j][i][uy];
            u_con[3] = x[j][i][uz];

            du_con_dt[0] = dgamma_dt;
            du_con_dt[1] = dx_dt[j][i][ux];
            du_con_dt[2] = dx_dt[j][i][uy];
            du_con_dt[3] = dx_dt[j][i][uz];

            u_cov[0] = -gamma;
            u_cov[1] = x[j][i][ux];
            u_cov[2] = x[j][i][uy];
            u_cov[3] = x[j][i][uz];

            du_cov_dt[0] = -dgamma_dt;
            du_cov_dt[1] = dx_dt[j][i][ux];
            du_cov_dt[2] = dx_dt[j][i][uy];
            du_cov_dt[3] = dx_dt[j][i][uz];

            b_con[0] = x[j][i][bx]*u_con[1] + x[j][i][by]*u_con[2] + x[j][i][bz]*u_con[3];
            b_con[1] = (x[j][i][bx] + b_con[0]*u_con[1])/u_con[0];
            b_con[2] = (x[j][i][by] + b_con[0]*u_con[2])/u_con[0];
            b_con[3] = (x[j][i][bz] + b_con[0]*u_con[3])/u_con[0];

            db_con_dt[0] = dx_dt[j][i][bx]*u_con[1] + dx_dt[j][i][by]*u_con[2] + 
                           dx_dt[j][i][bz]*u_con[3] + x[j][i][bx]*du_con_dt[1] + 
                           x[j][i][by]*du_con_dt[2] + x[j][i][bz]*du_con_dt[3];

            db_con_dt[1] = -(x[j][i][bx] + b_con[0]*u_con[1])*du_con_dt[0]/
                            (u_con[0]*u_con[0]) + 
                            (dx_dt[j][i][bx] + b_con[0]*du_con_dt[1] + 
                             db_con_dt[0]*u_con[1])/u_con[0];
        
            db_con_dt[2] = -(x[j][i][by] + b_con[0]*u_con[2])*du_con_dt[0]/
                            (u_con[0]*u_con[0]) + 
                            (dx_dt[j][i][by] + b_con[0]*du_con_dt[2] + 
                             db_con_dt[0]*u_con[2])/u_con[0];

            db_con_dt[3] = -(x[j][i][bz] + b_con[0]*u_con[3])*du_con_dt[0]/
                            (u_con[0]*u_con[0]) + 
                            (dx_dt[j][i][bz] + b_con[0]*du_con_dt[3] + 
                             db_con_dt[0]*u_con[3])/u_con[0];

            b_cov[0] = -b_con[0];
            b_cov[1] = b_con[1];
            b_cov[2] = b_con[2];
            b_cov[3] = b_con[3];

            db_cov_dt[0] = -db_con_dt[0];
            db_cov_dt[1] = db_con_dt[1];
            db_cov_dt[2] = db_con_dt[2];
            db_cov_dt[3] = db_con_dt[3];

            b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
                     b_con[3]*b_con[3];

            db_sqr_dt = 2.*(-b_con[0]*db_con_dt[0] + b_con[1]*db_con_dt[1] + 
                            b_con[2]*db_con_dt[2] + b_con[3]*db_con_dt[3]);

            dU_dt[j][i][rho] = dx_dt[j][i][rho]*u_con[dir] + x[j][i][rho]*du_con_dt[dir];

            dU_dt[j][i][bx] = dx_dt[j][i][bx];

            dU_dt[j][i][by] = dx_dt[j][i][by];

            dU_dt[j][i][bz] = dx_dt[j][i][bz];

//        dU_dt[j][i][bx] = db_con_dt[1]*u_con[dir] + b_con[1]*du_con_dt[dir] - 
//                       db_con_dt[dir]*u_con[1] - b_con[dir]*du_con_dt[1];
//
//        dU_dt[j][i][by] = db_con_dt[2]*u_con[dir] + b_con[2]*du_con_dt[dir] - 
//                       db_con_dt[dir]*u_con[2] - b_con[dir]*du_con_dt[2];
//
//        dU_dt[j][i][bz] = db_con_dt[3]*u_con[dir] + b_con[3]*du_con_dt[dir] - 
//                       db_con_dt[dir]*u_con[3] - b_con[dir]*du_con_dt[3];

            dU_dt[j][i][u] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
                           db_sqr_dt)*u_con[dir]*u_cov[0] + 
                          (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[0]
                          + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[0]
                          + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 0)
                          - db_con_dt[dir]*b_cov[0] - b_con[dir]*db_cov_dt[0];

            dU_dt[j][i][ux] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
                            db_sqr_dt)*u_con[dir]*u_cov[1] + 
                           (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[1]
                           + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[1]
                           + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 1)
                           - db_con_dt[dir]*b_cov[1] - b_con[dir]*db_cov_dt[1];
        
            dU_dt[j][i][uy] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
                            db_sqr_dt)*u_con[dir]*u_cov[2] + 
                           (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[2]
                           + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[2]
                           + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 2)
                           - db_con_dt[dir]*b_cov[2] - b_con[dir]*db_cov_dt[2];

            dU_dt[j][i][uz] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
                            db_sqr_dt)*u_con[dir]*u_cov[3] + 
                           (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[3]
                           + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[3]
                           + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 3)
                           - db_con_dt[dir]*b_cov[3] - b_con[dir]*db_cov_dt[3];

    }
    return(0);
}

PetscErrorCode ComputedU_dt3D(PetscScalar ****dU_dt, 
                              PetscScalar ****x, PetscScalar ****dx_dt,
                              AppCtx *ctx)
{
    PetscInt i, j, k, dir=0;
    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;
    PetscScalar dgamma_dt, du_cov_dt[4], du_con_dt[4];
    PetscScalar db_cov_dt[4], db_con_dt[4], db_sqr_dt, dP_dt;

    for (k=ctx->zstart; k<ctx->zstart + ctx->zsize; k++)
        for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
            for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++) {
                P = (4./3 - 1)*x[k][j][i][u];
                dP_dt = (4./3 - 1)*dx_dt[k][j][i][u];

                gamma = PetscSqrtScalar(1. + x[k][j][i][ux]*x[k][j][i][ux] + x[k][j][i][uy]*x[k][j][i][uy] +
                                             x[k][j][i][uz]*x[k][j][i][uz]);
                dgamma_dt = (x[k][j][i][ux]*dx_dt[k][j][i][ux] + x[k][j][i][uy]*dx_dt[k][j][i][uy] +
                            x[k][j][i][uz]*dx_dt[k][j][i][uz])/gamma;

                u_con[0] = gamma;
                u_con[1] = x[k][j][i][ux];
                u_con[2] = x[k][j][i][uy];
                u_con[3] = x[k][j][i][uz];

                du_con_dt[0] = dgamma_dt;
                du_con_dt[1] = dx_dt[k][j][i][ux];
                du_con_dt[2] = dx_dt[k][j][i][uy];
                du_con_dt[3] = dx_dt[k][j][i][uz];

                u_cov[0] = -gamma;
                u_cov[1] = x[k][j][i][ux];
                u_cov[2] = x[k][j][i][uy];
                u_cov[3] = x[k][j][i][uz];

                du_cov_dt[0] = -dgamma_dt;
                du_cov_dt[1] = dx_dt[k][j][i][ux];
                du_cov_dt[2] = dx_dt[k][j][i][uy];
                du_cov_dt[3] = dx_dt[k][j][i][uz];

                b_con[0] = x[k][j][i][bx]*u_con[1] + x[k][j][i][by]*u_con[2] + x[k][j][i][bz]*u_con[3];
                b_con[1] = (x[k][j][i][bx] + b_con[0]*u_con[1])/u_con[0];
                b_con[2] = (x[k][j][i][by] + b_con[0]*u_con[2])/u_con[0];
                b_con[3] = (x[k][j][i][bz] + b_con[0]*u_con[3])/u_con[0];

                db_con_dt[0] = dx_dt[k][j][i][bx]*u_con[1] + dx_dt[k][j][i][by]*u_con[2] + 
                               dx_dt[k][j][i][bz]*u_con[3] + x[k][j][i][bx]*du_con_dt[1] + 
                               x[k][j][i][by]*du_con_dt[2] + x[k][j][i][bz]*du_con_dt[3];

                db_con_dt[1] = -(x[k][j][i][bx] + b_con[0]*u_con[1])*du_con_dt[0]/
                                (u_con[0]*u_con[0]) + 
                                (dx_dt[k][j][i][bx] + b_con[0]*du_con_dt[1] + 
                                 db_con_dt[0]*u_con[1])/u_con[0];
        
                db_con_dt[2] = -(x[k][j][i][by] + b_con[0]*u_con[2])*du_con_dt[0]/
                                (u_con[0]*u_con[0]) + 
                                (dx_dt[k][j][i][by] + b_con[0]*du_con_dt[2] + 
                                 db_con_dt[0]*u_con[2])/u_con[0];

                db_con_dt[3] = -(x[k][j][i][bz] + b_con[0]*u_con[3])*du_con_dt[0]/
                                (u_con[0]*u_con[0]) + 
                                (dx_dt[k][j][i][bz] + b_con[0]*du_con_dt[3] + 
                                 db_con_dt[0]*u_con[3])/u_con[0];

                b_cov[0] = -b_con[0];
                b_cov[1] = b_con[1];
                b_cov[2] = b_con[2];
                b_cov[3] = b_con[3];

                db_cov_dt[0] = -db_con_dt[0];
                db_cov_dt[1] = db_con_dt[1];
                db_cov_dt[2] = db_con_dt[2];
                db_cov_dt[3] = db_con_dt[3];

                b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
                         b_con[3]*b_con[3];

                db_sqr_dt = 2.*(-b_con[0]*db_con_dt[0] + b_con[1]*db_con_dt[1] + 
                                b_con[2]*db_con_dt[2] + b_con[3]*db_con_dt[3]);

                dU_dt[k][j][i][rho] = dx_dt[k][j][i][rho]*u_con[dir] + x[k][j][i][rho]*du_con_dt[dir];

                dU_dt[k][j][i][bx] = dx_dt[k][j][i][bx];

                dU_dt[k][j][i][by] = dx_dt[k][j][i][by];

                dU_dt[k][j][i][bz] = dx_dt[k][j][i][bz];

                dU_dt[k][j][i][u] = (dP_dt + dx_dt[k][j][i][rho] + dx_dt[k][j][i][u] +
                                  db_sqr_dt)*u_con[dir]*u_cov[0] + 
                                 (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[0]
                                 + (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[0]
                                 + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 0)
                                 - db_con_dt[dir]*b_cov[0] - b_con[dir]*db_cov_dt[0];

                dU_dt[k][j][i][ux] = (dP_dt + dx_dt[k][j][i][rho] + dx_dt[k][j][i][u] +
                                   db_sqr_dt)*u_con[dir]*u_cov[1] + 
                                  (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[1]
                                  + (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[1]
                                  + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 1)
                                  - db_con_dt[dir]*b_cov[1] - b_con[dir]*db_cov_dt[1];
        
                dU_dt[k][j][i][uy] = (dP_dt + dx_dt[k][j][i][rho] + dx_dt[k][j][i][u] +
                                   db_sqr_dt)*u_con[dir]*u_cov[2] + 
                                  (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[2]
                                  + (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[2]
                                  + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 2)
                                  - db_con_dt[dir]*b_cov[2] - b_con[dir]*db_cov_dt[2];

                dU_dt[k][j][i][uz] = (dP_dt + dx_dt[k][j][i][rho] + dx_dt[k][j][i][u] +
                                   db_sqr_dt)*u_con[dir]*u_cov[3] + 
                                  (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[3]
                                  + (P + x[k][j][i][rho] + x[k][j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[3]
                                  + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 3)
                                  - db_con_dt[dir]*b_cov[3] - b_con[dir]*db_cov_dt[3];

    }
    return(0);
}

PetscErrorCode FluxCT2D(PetscScalar ***flux_x,
                        PetscScalar ***flux_y, AppCtx *ctx)
{
    if (ctx->dim!=2) 
        SETERRQ1(PETSC_COMM_SELF, 1, "2D FluxCT called for %d dim\n", ctx->dim);

    PetscInt i, j;

    for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++)
        for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++) {
            flux_x[j][i][bx] = 0.0;
            
            flux_y[j][i][by] = 0.0;

            flux_y[j][i][bx] = 0.125*(2.*flux_y[j][i][bx] +
                                      flux_y[j][i+1][bx] + 
                                      flux_y[j][i-1][bx] -
                                      flux_x[j][i][by] -
                                      flux_x[j][i+1][by] -
                                      flux_x[j-1][i][by] -
                                      flux_x[j-1][i+1][by]);

            flux_x[j][i][by] = 0.125*(2.*flux_x[j][i][by] +
                                      flux_x[j+1][i][by] +
                                      flux_x[j-1][i][by] -
                                      flux_y[j][i][bx] -
                                      flux_y[j+1][i][bx] -
                                      flux_y[j][i-1][bx] -
                                      flux_y[j+1][i-1][bx]);
        }

    return(0);

}

PetscErrorCode RiemannSolver1D(PetscScalar **flux_x, 
                               PetscScalar **flux_lx, 
                               PetscScalar **flux_rx,
                               PetscScalar **U_lx,
                               PetscScalar **U_rx, AppCtx *ctx)
{
    if (ctx->dim!=1) 
        SETERRQ1(PETSC_COMM_SELF, 1, "1D RiemannSolver called for %d dim\n", ctx->dim);

    PetscInt i, var;

    for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++)
        for (var=0; var<dof; var++)
            flux_x[i][var] = 0.5*(flux_lx[i][var] + flux_rx[i-1][var] -
                                  (U_lx[i][var] - U_rx[i-1][var]));

    return(0);

}

PetscErrorCode RiemannSolver2D(PetscScalar ***flux_x, 
                               PetscScalar ***flux_y,
                               PetscScalar ***flux_lx, 
                               PetscScalar ***flux_rx,
                               PetscScalar ***flux_ly, 
                               PetscScalar ***flux_ry,
                               PetscScalar ***U_lx, 
                               PetscScalar ***U_rx,
                               PetscScalar ***U_ly, 
                               PetscScalar ***U_ry, AppCtx *ctx)
{
    if (ctx->dim!=2) 
        SETERRQ1(PETSC_COMM_SELF, 1, "2D RiemannSolver called for %d dim\n", ctx->dim);

    PetscInt i, j, var;

    for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++)
        for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++)
            for (var=0; var<dof; var++) {
                flux_x[j][i][var] = 0.5*(flux_lx[j][i][var] + flux_rx[j][i-1][var] -
                                        (U_lx[j][i][var] - U_rx[j][i-1][var]));
                flux_y[j][i][var] = 0.5*(flux_ly[j][i][var] + flux_ry[j-1][i][var] -
                                        (U_ly[j][i][var] - U_ry[j-1][i][var]));
            }

    return(0);

}

PetscErrorCode RiemannSolver3D(PetscScalar ****flux_x, 
                               PetscScalar ****flux_y,
                               PetscScalar ****flux_z,
                               PetscScalar ****flux_lx, 
                               PetscScalar ****flux_rx,
                               PetscScalar ****flux_ly, 
                               PetscScalar ****flux_ry,
                               PetscScalar ****flux_lz, 
                               PetscScalar ****flux_rz,
                               PetscScalar ****U_lx, 
                               PetscScalar ****U_rx,
                               PetscScalar ****U_ly, 
                               PetscScalar ****U_ry,
                               PetscScalar ****U_lz, 
                               PetscScalar ****U_rz, AppCtx *ctx)
{
    if (ctx->dim!=3) 
        SETERRQ1(PETSC_COMM_SELF, 1, "3D RiemannSolver called for %d dim\n", ctx->dim);

    PetscInt i, j, k, var;

    for (k=ctx->zstart-1; k<ctx->zstart + ctx->zsize+1; k++)
        for (j=ctx->ystart-1; j<ctx->ystart + ctx->ysize+1; j++)
            for (i=ctx->xstart-1; i<ctx->xstart + ctx->xsize+1; i++)
                for (var=0; var<dof; var++) {

                    flux_x[k][j][i][var] = 0.5*(flux_lx[k][j][i][var] + 
                                                flux_rx[k][j][i-1][var] -
                                               (U_lx[k][j][i][var] - 
                                                U_rx[k][j][i-1][var]));

                    flux_y[k][j][i][var] = 0.5*(flux_ly[k][j][i][var] +
                                                flux_ry[k][j-1][i][var] -
                                               (U_ly[k][j][i][var] - 
                                                U_ry[k][j-1][i][var]));

                    flux_z[k][j][i][var] = 0.5*(flux_lz[k][j][i][var] +
                                                flux_rz[k-1][j][i][var] -
                                               (U_lz[k][j][i][var] - 
                                                U_rz[k-1][j][i][var]));
                }

    return(0);
}

PetscErrorCode ResFunction(TS ts,
                           PetscScalar t, 
                           Vec Xvec, Vec dX_dt, Vec F,
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
            PetscScalar **x, **f, **dx_dt;
            PetscScalar **prim_lx, **prim_rx;
            PetscScalar **flux_lx, **flux_rx;
            PetscScalar **flux_x;
            PetscScalar **U_lx, **U_rx;

            DMDAVecGetArrayDOF(da, localX, &x);
            DMDAVecGetArrayDOF(da, F, &f);
            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);

            DMDAVecGetArrayDOF(da, ctx->prim_l[X], &prim_lx);
            DMDAVecGetArrayDOF(da, ctx->prim_r[X], &prim_rx);

            DMDAVecGetArrayDOF(da, ctx->flux_l[X], &flux_lx);
            DMDAVecGetArrayDOF(da, ctx->flux_r[X], &flux_rx);

            DMDAVecGetArrayDOF(da, ctx->flux[X], &flux_x);

            DMDAVecGetArrayDOF(da, ctx->U_l[X], &U_lx);
            DMDAVecGetArrayDOF(da, ctx->U_r[X], &U_rx);

            ComputedU_dt1D(f, x, dx_dt, ctx);

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
                    f[i][var] = f[i][var] + 
                                (flux_x[i+1][var] - flux_x[i][var])/ctx->dx;

            DMDAVecRestoreArrayDOF(da, localX, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);
            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);

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
            PetscScalar ***x, ***f, ***dx_dt;
            PetscScalar ***prim_lx, ***prim_rx;
            PetscScalar ***prim_ly, ***prim_ry;
            PetscScalar ***flux_lx, ***flux_rx;
            PetscScalar ***flux_ly, ***flux_ry;
            PetscScalar ***flux_x, ***flux_y;
            PetscScalar ***U_lx, ***U_rx;
            PetscScalar ***U_ly, ***U_ry;

            DMDAVecGetArrayDOF(da, localX, &x);
            DMDAVecGetArrayDOF(da, F, &f);
            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);

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

            ComputedU_dt2D(f, x, dx_dt, ctx);

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

            FluxCT2D(flux_x, flux_y, ctx);

            for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
                for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++)
                    for (var=0; var<dof; var++) {
                        f[j][i][var] =  f[j][i][var] + 
                                        (flux_x[j][i+1][var] - 
                                        flux_x[j][i][var])/ctx->dx +
                                        (flux_y[j+1][i][var] -
                                        flux_y[j][i][var])/ctx->dy;
                    }


            DMDAVecRestoreArrayDOF(da, localX, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);
            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);

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
            PetscScalar ****x, ****f, ****dx_dt;
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
            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);

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

            ComputedU_dt3D(f, x, dx_dt, ctx);
            DMGlobalToLocalEnd(da, Xvec, INSERT_VALUES, localX);

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
                            f[k][j][i][var] = f[k][j][i][var] + 
                                            (flux_x[k][j][i+1][var] - 
                                            flux_x[k][j][i][var])/ctx->dx +
                                           
                                            (flux_y[k][j+1][i][var] -
                                            flux_y[k][j][i][var])/ctx->dy +

                                            (flux_z[k+1][j][i][var] -
                                            flux_z[k][j][i][var])/ctx->dz;

                        }

            DMDAVecRestoreArrayDOF(da, localX, &x);
            DMDAVecRestoreArrayDOF(da, F, &f);
            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);

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
