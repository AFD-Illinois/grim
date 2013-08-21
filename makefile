ALL: grim

PETSC_DIR=/home/mc/Downloads/petsc_optimized

#CFLAGS = -g
#CFLAGS = -ftree-vectorizer-verbose=2

FFLAGS =

#CPPFLAGS = -g
#CPPFLAGS = -ftree-vectorizer-verbose=2

FPPFLAGS =

include ${PETSC_DIR}/conf/variables

include ${PETSC_DIR}/conf/rules

grim: grim.o
	-${CLINKER} -o grim grim.o ${PETSC_LIB}
	${RM} grim.o
