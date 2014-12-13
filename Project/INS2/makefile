ALL: ibm2d

SRC = ${wildcard *.cpp}
OBJ = ${SRC: .cpp=.o}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I.
PCC_FLAGS += -std=c++0x -Wall -Wextra -pedantic

NS2d.o: main.cpp NavierStokesSolver.cpp fluidsimulation.cpp
	${PCC} -g -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} $^

ibm2d: ${OBJ}
	${CLINKER} -g $^ -o $@ ${PETSC_SYS_LIB}
	
