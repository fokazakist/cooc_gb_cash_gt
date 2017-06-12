CXX = g++
VERSION = 0.01
NAME = gspan
CXXFLAGS = -O2 -Wall -Wno-deprecated
LDFLAGS =  -lm -ldl #remove "-static" for Mac,in linux need to compile?
TARGETS = lpboost
OBJ = gspan.o cashtree.o 

all:  lpboost 

gspan: gspan.o
	${CXX} ${LDFLAGS} -o ${TARGETS} ${OBJ}  ${LDFLAGS}

lpboost: lpboost.o ${OBJ}
	${CXX} ${CXXFLAGS} -o ${TARGETS} ${OBJ} lpboost.o ${LDFLAGS} -I../glpk-4.8/include/ ../glpk-4.8/src/libglpk.a

clean:
	rm -f ${OBJ} lpboost.o ${TARGETS} *~

