CXX     = g++
MPICXX  = mpic++
LIB     = -lmpi -lpthread -lgfortran -llapack -lblas
FLAG    = -Wall -I$(INCLUDE) -O2 -fpermissive
DEF_OPENMP = -fopenmp -DCAPLET_OPENMP
DEF_MPI = -DCAPLET_MPI

SRC        = src
OBJ_MPI    = obj/mpi
OBJ_OPENMP = obj/openmp
BIN        = bin
INCLUDE    = include
EXAMPLE    = example

CAPLET_MPI_OBJ = \
	$(OBJ_MPI)/caplet.o \
	$(OBJ_MPI)/caplet_elem.o \
	$(OBJ_MPI)/caplet_int.o \
	$(OBJ_MPI)/caplet_widgets.o \
	$(OBJ_MPI)/main.o \

CAPLET_OPENMP_OBJ = \
	$(OBJ_OPENMP)/caplet.o \
	$(OBJ_OPENMP)/caplet_elem.o \
	$(OBJ_OPENMP)/caplet_int.o \
	$(OBJ_OPENMP)/caplet_widgets.o \
	$(OBJ_OPENMP)/main.o \

CAPLET_INCLUDE = \
	$(INCLUDE)/caplet.h \
	$(INCLUDE)/caplet_blas.h \
	$(INCLUDE)/caplet_const.h \
	$(INCLUDE)/caplet_debug.h \
	$(INCLUDE)/caplet_elem.h \
	$(INCLUDE)/caplet_gauss.h \
	$(INCLUDE)/caplet_int.h \
	$(INCLUDE)/caplet_parameter.h \
	$(INCLUDE)/caplet_widgets.h

all: capletMPI capletOpenMP 

mpi: capletMPI 

openmp: capletOpenMP 

capletMPI: $(CAPLET_MPI_OBJ)
	$(MPICXX) $(FLAG) $(DEF_MPI) -o $(BIN)/$@ $^ $(LIB)

$(OBJ_MPI)/%.o: $(SRC)/%.cpp $(CAPLET_INCLUDE)
	$(MPICXX) $(FLAG) $(DEF_MPI) -c $< -o $@

capletOpenMP: $(CAPLET_OPENMP_OBJ)
	$(MPICXX) $(FLAG) $(DEF_OPENMP) -o $(BIN)/$@ $^ $(LIB) 

$(OBJ_OPENMP)/%.o: $(SRC)/%.cpp $(CAPLET_INCLUDE)
	$(MPICXX) $(FLAG) $(DEF_OPENMP) -c $< -o $@

.phony: clean
clean:
	rm -rf $(OBJ_MPI)/* $(OBJ_OPENMP)/* 
