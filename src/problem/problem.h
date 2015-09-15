#ifndef GRIM_PROBLEM_H_
#define GRIM_PROBLEM_H_

#include <petsc.h>
#include "../inputs.h"
#include "../timestepper/timestepper.h"
#include "../gridzone/gridzone.h"
#include "../geometry/geometry.h"
#include "../boundary/boundary.h"
#include "../physics/physics.h"
#include PROBLEM_DATA /* Includes the problem struct defined in the specific
                         problem folder */

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);

void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE]);

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
                                     const int X1Start, const int X2Start,
                                     const int X1Size, const int X2Size,
                                     REAL primTile[ARRAY_ARGS TILE_SIZE]);

void applyAdditionalProblemSpecificBCs
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL primTile[ARRAY_ARGS TILE_SIZE]
);

void applyProblemSpecificFluxFilter
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
);

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

void writeProblemSpecificData(PetscViewer parametersViewer,
    const struct problemData problemSpecificData[ARRAY_ARGS 1]);

#if (CONDUCTION)
  void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1]);
#endif
#if (VISCOSITY)
  void setViscosityParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1]);
#endif

#define WRITE_PARAM_INT(NAME) \
        do { \
        int __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_INT, &__NAME); \
        } while(0)


#define WRITE_PARAM_DOUBLE(NAME) \
        do { \
        double __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_DOUBLE, &__NAME); \
        } while(0)


#define WRITE_PARAM_COMPLEX(NAME) \
        do { \
        complex __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_COMPLEX, &__NAME); \
        } while(0)


#define WRITE_PARAM_LONG(NAME) \
        do { \
        long __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_LONG, &__NAME); \
        } while(0)


#define WRITE_PARAM_SHORT(NAME) \
        do { \
        short __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_SHORT, &__NAME); \
        } while(0)


#define WRITE_PARAM_FLOAT(NAME) \
        do { \
        float __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_FLOAT, &__NAME); \
        } while(0)


#define WRITE_PARAM_CHAR(NAME) \
        do { \
        char __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_CHAR, &__NAME); \
        } while(0)


#define WRITE_PARAM_STRING(NAME) \
        do { \
        char *__NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_STRING, &__NAME); \
        } while(0)

#define WRITE_PARAM_INT(NAME) \
        do { \
        int __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_INT, &__NAME); \
        } while(0)


#define WRITE_PARAM_DOUBLE(NAME) \
        do { \
        double __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_DOUBLE, &__NAME); \
        } while(0)


#define WRITE_PARAM_COMPLEX(NAME) \
        do { \
        complex __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_COMPLEX, &__NAME); \
        } while(0)


#define WRITE_PARAM_LONG(NAME) \
        do { \
        long __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_LONG, &__NAME); \
        } while(0)


#define WRITE_PARAM_SHORT(NAME) \
        do { \
        short __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_SHORT, &__NAME); \
        } while(0)


#define WRITE_PARAM_FLOAT(NAME) \
        do { \
        float __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_FLOAT, &__NAME); \
        } while(0)


#define WRITE_PARAM_CHAR(NAME) \
        do { \
        char __NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_CHAR, &__NAME); \
        } while(0)


#define WRITE_PARAM_STRING(NAME) \
        do { \
        char *__NAME = NAME; \
        const char* __NAME_GROUP; \
        PetscViewerHDF5GetGroup(parametersViewer, &__NAME_GROUP); \
        PetscViewerHDF5WriteAttribute(parametersViewer, \
        __NAME_GROUP, #NAME, PETSC_STRING, &__NAME); \
        } while(0)

#endif /* GRIM_PROBLEM_H_ */
