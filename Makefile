#
# A code for adaptively-smoothed density, velocity, and clumping fields
#

EXEC      = clumping_tree

MPICXX   ?= mpicxx
OPTIMIZE ?= -O2
OPENMP    = -openmp
OPTIONS  := $(OPTIMIZE) $(OPENMP)

OPTIONS  += -DHALO_EXCISE # Excise halo using given halo catalog
# OPTIONS  += FOFON  # Do also Friends-of-friends halo finding/excising
# OPTIONS  += -DBGQ  # For big endian machine (file_reader2.cpp)
OPTIONS  += -DPID_FLAG # Halo file format with particle IDs


OBJS     := clumping_tree_main.o 
OBJS     += option.o
OBJS     += basic_types.o mpi_interface.o
OBJS     += logger.o
OBJS     += file_reader2.o exchange_buffer.o
OBJS     += kdtree_balanced.o
OBJS     += nbr_search.o fof.o coarse_mesh.o open.o
OBJS     += halo_file2.o
OBJS     += global_mesh.o

CXXFLAGS  = $(OPTIONS)

LIBS      = -lm

all: $(EXEC)
clumping_tree: $(OBJS) performance_items.h
	$(MPICXX) $(OPENMP) $(OBJS) $(LIBS) -o $@

performance_items.h: performance_items.sh performance_items.txt 
	sh $< > $@

%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<




.PHONY : clean test
clean:
	rm -f $(OBJS) $(EXEC)






