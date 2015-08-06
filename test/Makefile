ROOT = ..

include ${ROOT}/Makefile

TESTNAME=testbiharmonic
TESTFILE=testBiharmonic.cpp

test:
	rm ${TESTNAME}; \
	cd ${ROOT}/src; \
        make objs; \
	cd ${ROOT}/test; \
	make ${TESTNAME};

${TESTNAME}:
	${CC} -o ${TESTNAME} ${TESTFILE} ${PS_OBJ} \
        ${HDF5_LIB} ${CFLAGS} ${PS_LIB} ${PS_INC}
