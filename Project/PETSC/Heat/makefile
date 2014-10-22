ALL:heat_petsc

CFLGAS =
FFLAGS =
CPPFLAGS =
FPPFLAGS =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
#include ${PETSC_DIR}/bmake/common/base

heat_petsc:heat_petsc.o chkopts
	-${CLINKER} -o heat_petsc heat_petsc.o ${PETSC_SNES_LIB}
	${RM} heat_petsc.o

