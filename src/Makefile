ROOT = ..

include $(ROOT)/Makefile

objs: Field.o \
      Force.o \
      Utils.o \
      QGEqn.o \
      libbase.a

QGEqn.o: QGEqn.cpp QGEqn.h Force.o Utils.o Field.o
	${CC} -c QGEqn.cpp ${CFLAGS}
Force.o: Force.cpp Force.h Field.o
	${CC} -c Force.cpp ${CFLAGS}
Utils.o: Utils.cpp Utils.h Field.o
	${CC} -c Utils.cpp ${CFLAGS}
Field.o: Field.cpp Field.h
	${CC} -c Field.cpp ${CFLAGS}

libbase.a:
	ar -cq libbase.a *.o
	mv libbase.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
