#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -g  
OPTFLAGS   = -O3  
#OPTFLAGS   = 

LINKER      = g++
LDFLAGS   = -lgsl -lgslcblas -lm -lboost_program_options-mt

#------ Intel icpc  --------------------------------------------------

#CC         = icpc
#CPPFLAGS     = -g -Wall
#OPTFLAGS   = -O3 -ipo -fnsplit -ip -ftz -funroll-loops
#OPTFLAGS   = -static
#OPTFLAGS   = -fast / -mp / 
#OPTFLAGS   = 

#LINKER      = icpc
#LDFLAGS   = -lm -lscs -lmkl 
#LDFLAGS   = -i-dynamic/static
#LDFLAGS     = 

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

#SOURCES     = mfa_solve.cc ode.cc ode_io_utils.cc mfa.cc 
#OBJECTS     = $(SOURCES:.cc=.o)

SOURCES     = solve_ode.cc 
OBJECTS     = $(SOURCES:.cc=.o)

EXEC        = solve_ode.x
#--------------------------------------------------------------------

all:  $(EXEC)

$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

%.o : %.cc
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) $(SOURCES)

clean:
	@echo ... cleaning $(OBJECTS)
	@rm -f $(OBJECTS)
	@echo ... done

cleaner:
	@echo '... cleaning also *~, *.x, gmon.out and '$(EXEC)
	@rm -f $(OBJECTS) *~ *.x gmon.out $(EXEC)
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

solve_ode.o: file_io.hh Sleep.hh pa_macros.hh \
             ode_solver.hh ode_io_utils.hh convergence.hh \
             Diminish.hh
