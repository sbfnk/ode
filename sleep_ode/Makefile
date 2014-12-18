#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -g -O3 -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS   = -lgsl -lgslcblas -lm -lboost_program_options-mt

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = solve_ode.cc 
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = solve_ode.x
#--------------------------------------------------------------------

all:  bin/$(EXEC)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

out/%.o : src/%.cc include/%.hh
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) $(SOURCES)

clean:
	@echo ... cleaning $(OBJECTS) $(EXEC)
	@rm -f bin/* out/*
	@echo ... done

cleaner:
	@echo '... cleaning also *~, *.x and gmon.out
	@rm -f bin/* out/* *~ *.x gmon.out $(EXEC)
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

solve_ode.o: file_io.hh Sleep.hh ode_solver.hh ode_io_utils.hh 
	     convergence.hh 
