all:example
	#if there are additional programs, their name would be placed after example
	include ${PETSC_DIR}
	#library required duirng compilation
	src=example.c utilities.c parameters.c guess.c equations.c jacobian.c
	obj=example.o utilities.o parameters.o guess.o equations.o jacobian.o
	# "make" will scan the source files for changes since last complie and recompile if necessary
	
	#linking 
example:${obj} chkopts
	${CLINKER} -o $@ $< ${PETSC_LIB}
	${RM} example.o
	#${PETSC_LIB} link object files with all of PETSc libs


