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

# Particle file format
# OPTIONS  += -DBGQ  # Running this on a big endian machine (file_reader2.cpp)
# OPTIONS += -DSWAPENDIAN # Running this code on different endian wrt file
OPTIONS  += -DPID_FLAG # Halo file format with particle IDs

# Halo removal
# OPTIONS  += -DM200 # Use M200/r200 instead of Mvir/rvir (Bryan and Norman)

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

# % g++ -MM *.cpp
$(OBJS): Makefile

basic_types.o: basic_types.cpp basic_types.h
clumping_tree_main.o: clumping_tree_main.cpp option.h logger.h \
  performance_items.h basic_types.h mpi_interface.h file_reader2.h \
  exchange_buffer.h kdtree_balanced.h nbr_search.h k_neighbors.h fof.h \
  coarse_mesh.h halo_file.h open.h
coarse_mesh.o: coarse_mesh.cpp basic_types.h coarse_mesh.h open.h
exchange_buffer.o: exchange_buffer.cpp exchange_buffer.h basic_types.h \
  mpi_interface.h
file_reader.o: file_reader.cpp basic_types.h file_reader.h \
  mpi_interface.h
file_reader2.o: file_reader2.cpp basic_types.h file_reader.h \
  mpi_interface.h
fof.o: fof.cpp fof.h basic_types.h kdtree_balanced.h
global_mesh.o: global_mesh.cpp basic_types.h mpi_interface.h \
  global_mesh.h open.h
halo_file.o: halo_file.cpp kdtree_balanced.h basic_types.h halo_file.h \
  mpi_interface.h exchange_buffer_template.h nbr_search.h k_neighbors.h
halo_file2.o: halo_file2.cpp kdtree_balanced.h basic_types.h halo_file.h \
  mpi_interface.h exchange_buffer_template.h nbr_search.h k_neighbors.h
histogram.o: histogram.cpp histogram.h
kdtree_balanced.o: kdtree_balanced.cpp kdtree_balanced.h basic_types.h
logger.o: logger.cpp logger.h performance_items.h
mpi_interface.o: mpi_interface.cpp mpi_interface.h
nbr_search.o: nbr_search.cpp nbr_search.h basic_types.h kdtree_balanced.h \
  k_neighbors.h
open.o: open.cpp open.h
option.o: option.cpp option.h



.PHONY : clean test
clean:
	rm -f $(OBJS) $(EXEC)
