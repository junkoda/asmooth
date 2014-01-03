#
# A code for adaptively-smoothed density, velocity, and clumping fields
#

EXEC      = clumping_tree #clumping_tree_async density_pdf fof_only

MPICXX   ?= mpicxx
OPTIMIZE ?=  -O2
OPENMP    = -openmp
OPTIONS   := $(OPTIMIZE) -I$(MYLIB) $(OPENMP)

OPTIONS  += HALO_EXCISE # Excise halo using given halo catalog
# OPTIONS  += FOFON  # Do also Friends-of-friends halo finding/excising

OBJS     := option.o
OBJS     += basic_types.o mpi_interface.o
OBJS     += logger.o
OBJS     += file_reader.o exchange_buffer.o
OBJS     += kdtree_balanced.o
OBJS     += nbr_search.o fof.o coarse_mesh.o open.o
OBJS     += halo_file.o
OBJS     += global_mesh.o

OBJS2     = fof_only_main.o option.o basic_types.o mpi_interface.o logger.o file_reader.o exchange_buffer.o kdtree_balanced.o nbr_search.o fof.o coarse_mesh.o open.o

INCL      = architecture_params.h basic_types.h logger.h mpi_interface.h performance_items.h file_reader.h exchange_buffer.h kdtree_balanced.h k_neighbors.h nbr_search.h fof.h coarse_mesh.h open.h


CXXFLAGS  = $(OPTIONS)

LIBS      = -lm -lstdc++ $(PTHREAD)

all: $(EXEC)
clumping_tree: clumping_tree_main.o $(OBJS) 
	$(MPICXX) $(OPENMP) $< $(OBJS) $(LIBS) -o $@

clumping_tree_async: clumping_tree_main_asynchronous.o $(OBJS)
	$(MPICXX) $(OPENMP) $< $(OBJS) $(LIBS) -o $@

density_pdf: density_pdf_main.o histogram.o $(OBJS)
	$(MPICXX) $(OPENMP) $< histogram.o $(OBJS) $(LIBS) -o $@

performance_items.h: performance_items.sh performance_items.txt 
	sh $< > $@

fof_only: $(OBJS2)
	$(MPICXX) $(OPENMP) $(OBJS2) $(LIBS) -o $@

%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<


$(OBJS): $(INCL) 
clumping_tree_main.o: $(INCL)
clumping_tree_main_asynchronous.o: $(INCL)
density_pdf_main.o: $(INCL)
$(OBJS2): $(INCL)

logger.o: logger.h performance_items.h



.PHONY : clean test
clean:
	rm -f $(OBJS) clumping_tree_main.o clumping_tree_main_asynchronous.o density_pdf_main.o histogram.o $(EXEC)

#	cd test && $(MAKE) $@
#test:
#	cd $@ && $(MAKE) test








