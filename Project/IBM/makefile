ALL: ibm2d

#SRC = ${wildcard *.cpp}
#OBJ = ${SRC: .cpp=.o}

SRC = main.cpp TcSolver.cpp INS2/NavierStokesSolver.cpp
OBJ = ${SRC: .cpp = .o}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_CC_INCLUDES += -I. -I./INS2
# PCC_FLAGS += -std=c++0x -Wall -Wextra -pedantic
PCC_FLAGS += -std=gnu++0x -Wall -Wextra -pedantic

ibm2d.o: ${SRC} 
	${PCC} -g -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} $^

ibm2d: ${OBJ}
	${CLINKER} -g $^ -o $@ ${PETSC_SYS_LIB}
	
