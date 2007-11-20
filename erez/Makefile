#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS     = -lboost_program_options -lgsl 

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = solve_ode.cc
OBJECTS     = $(patsubst %.cc,obj/%.o,${SOURCES})

EXEC        = solve_ode.x
#--------------------------------------------------------------------

all:  bin/$(EXEC)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking: 
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

obj/%.o : src/%.cc include/%.hh
#out/%.o : src/%.cc 
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) src/*.cc

htmldoc:
	doxygen

pdfdoc: htmldoc
	$(MAKE) -C doc/latex pdf

clean:
	@echo ... cleaning $(OBJECTS) $(EXEC)
	@rm -f bin/* obj/*
	@echo ... done

cleaner:
	@echo '... cleaning also *~'
	@rm -f bin/* obj/* *~ .#*
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

obj/solve_ode.o: include/ode_solver.hh include/InfoSIRSmf.hh \
                 include/InfoSIRSpa.hh include/ode_io_utils.hh \
		 include/convergence.hh include/pa_macros.hh

