ROOT = ..

include ${ROOT}/Makefile

TOOLNAME=energySpectra
TOOLFILE=EnergySpectra.cpp

tool:
	rm ${TOOLNAME}; \
	cd ${ROOT}/src; \
        make objs; \
	cd ${ROOT}/tool; \
	make ${TOOLNAME};

${TOOLNAME}:
	${CC} -o ${TOOLNAME} ${TOOLFILE} ${PS_OBJ} \
        ${HDF5_LIB} ${CFLAGS} ${PS_LIB} ${PS_INC}
