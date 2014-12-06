# Dump folder
set(OUTPUT_DIR "\"/home/mc/data/mc/PhD/grim/build/\"")
set(RESTART "OFF")
set(RESTART_FILE "\"restartfile.h5\"")
set(DUMP_FILE_PREFIX "\"data\"")

# Time stepping options: EXPLICIT, IMEX or IMPLICIT
set(TIME_STEPPING "EXPLICIT")
set(DT "0.001")
set(DT_DUMP ".1")
set(START_TIME "0.")
set(FINAL_TIME "4.")
set(COURANT "0.5")

# Domain size. If the problem is 1D then N2 is ignored.
set(COMPUTE_DIM "2")
set(N1 "64")
set(N2 "64")

# The entire global domain is divided into tiles that are small enough to fit
# into the cache of the compute node or an accelerator. This technique
# optimizes cache usage and makes explicit use of the deep memory hierarchies
# which are prevalent in all kinds of processors and accelerators. 
# 
# Caution: 1) The tile sizes need to divide the global domain size in each
#             direction exactly!
#          2) We want a tile size that is as large as the cache. Too large a
#             size and the code can fail silently because local memory in the
#             OpenCL specification is a different amount for different 
#             machines.
#             NEED TO PUT A CHECK FOR THIS.
set(TILE_SIZE_X1 "8")
set(TILE_SIZE_X2 "8")

# Physics variables
set(ADIABATIC_INDEX "5./3")
set(CONDUCTION "OFF")

# Geometry
set(METRIC "MINKOWSKI")
set(EPS    "1e-5")

# Domain
set(X1_A  "-6.")
set(X1_B  "6.")
set(X2_A  "-6.")
set(X2_B  "6.")

# Boundary conditions
set(PHYSICAL_BOUNDARY_LEFT_EDGE   "OUTFLOW")
set(PHYSICAL_BOUNDARY_RIGHT_EDGE  "OUTFLOW")
set(PHYSICAL_BOUNDARY_TOP_EDGE    "OUTFLOW")
set(PHYSICAL_BOUNDARY_BOTTOM_EDGE "OUTFLOW")

# Reconstrution options
# MP5 or MONOTONIZED_CENTRAL
set(RECONSTRUCTION "MONOTONIZED_CENTRAL")

# Initial condition parameters 
# MEAN_MAGNETIC_FIELD -- Magnetic field in the x direction
# SMOOTHNESS -- Parameter for the shell edge
# Parameters taken from Komissarov, 1999
set(INITIAL_RADIUS      "1.")
set(MEAN_MAGNETIC_FIELD "0.")
set(RHO_INSIDE          "1e-2")
set(RHO_OUTSIDE         "1e-4")
set(PRESSURE_INSIDE     "1.")
set(PRESSURE_OUTSIDE    "3e-5")
set(SMOOTHNESS          "20.")

# Floor values
set(RHO_FLOOR "1e-7")
SET(UU_FLOOR "1e-9")

message("")
message("#########################")
message("# Configuration options #")
message("#########################")
message("")
message("Problem          : " ${PROBLEM})
message("Output dir       : " ${OUTPUT_DIR})
message("Time stepping    : " ${TIME_STEPPING})
message("Dimensions       : " ${COMPUTE_DIM})
message("Resolution       : " ${N1} " x " ${N2})
message("Tile size        : " ${TILE_SIZE_X1} " x " ${TILE_SIZE_X2})
message("Metric           : " ${METRIC})

message("")
message("##################################")
message("# Domain and boundary conditions #")
message("##################################")
message("")
message("              " ${PHYSICAL_BOUNDARY_TOP_EDGE})
message("    (" ${X1_A},${X2_B} ")+-----------+" "(" ${X1_B},${X2_B} ")" )
message(                "           |           |"       )
message(                "           |           |"       )
message("  " ${PHYSICAL_BOUNDARY_LEFT_EDGE} " |           | " ${PHYSICAL_BOUNDARY_RIGHT_EDGE})
message(                "           |           |"       )
message(                "           |           |"       )
message("    (" ${X1_A},${X2_A} ")+-----------+" "(" ${X1_B},${X2_A} ")" )
message("              " ${PHYSICAL_BOUNDARY_BOTTOM_EDGE})

message("")
message("###################")
message("# Physics options #")
message("###################")
message("")
message("Adiabatic index  : " ${ADIABATIC_INDEX})
message("Conduction       : " ${CONDUCTION})

message("")
message("##################")
message("# Reconstruction #")
message("##################")
message("")
message("Reconstruction  : " ${RECONSTRUCTION})

message("")
message("######################")
message("# Initial conditions #")
message("######################")
message("")
message("Initial radius       : " ${INITIAL_RADIUS})
message("Mean magnetic field  : " ${MEAN_MAGNETIC_FIELD})
message("Density inside       : " ${RHO_INSIDE})
message("Density outside      : " ${RHO_OUTSIDE})
message("Pressure inside      : " ${PRESSURE_INSIDE})
message("Pressure outside     : " ${PRESSURE_OUTSIDE})
message("Smoothness factor    : " ${SMOOTHNESS})

message("")
message("##########")
message("# Floors #")
message("##########")
message("")
message("Density floor            : " ${RHO_FLOOR})
message("Internal energy floor    : " ${UU_FLOOR})
