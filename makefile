ALL: grim

#PETSC_DIR=/home/mc/Downloads/petsc_icc
PETSC_DIR=/home/mc/Downloads/petsc_optimized

#CFLAGS = -O3 -xHOST -ipo -no-prec-div -xsse4.2 -msse4.2 -mkl -parallel -vec-report2 -align
CFLAGS = -ftree-vectorizer-verbose=2
#CFLAGS = -g

FFLAGS =

#CPPFLAGS = -O3 -xHOST -ipo -no-prec-div -xsse4.2 -msse4.2 -mkl -parallel -vec-report2 -align
CPPFLAGS = -ftree-vectorizer-verbose=2
#CPPFLAGS = -g

FPPFLAGS =

include ${PETSC_DIR}/conf/variables

include ${PETSC_DIR}/conf/rules

grim: grim.o
	-${CLINKER} -o grim grim.o ${PETSC_LIB}
	${RM} grim.o
