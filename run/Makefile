ROOT = ..

include ${ROOT}/Makefile

run:
	rm qgmodel; \
	cd ${ROOT}/src; \
        make objs; \
	cd ${ROOT}/run; \
	make qgmodel;

qgmodel:
	${CC} -o qgmodel qgModel_Forced.cpp ${PS_OBJ} \
        ${HDF5_LIB} ${CFLAGS} ${PS_LIB} ${PS_INC}

genpv:
	${CC} -o genpv GenPV.cpp ${PS_OBJ} \
        ${HDF5_LIB} ${CFLAGS} ${PS_LIB} ${PS_INC}
