#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -g  
#OPTFLAGS   = -O3 -ipo -tpp2 -fnsplit -ip -ftz 
OPTFLAGS   = 

LINKER      = g++
LDFLAGS   = -lgsl -lcblas -lm 

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

SOURCES     = ode_solve.cc ode.cc ode_io_utils.cc model_ode.cc 
OBJECTS     = $(SOURCES:.cc=.o)

EXEC        = out
#--------------------------------------------------------------------

all:  $(EXEC)

$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@.x
	@echo

%.o : %.cc
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) *.cc

clean:
	@echo ... cleaning $(OBJECTS)
	@rm -f $(OBJECTS)
	@echo ... done

cleaner:
	@echo '... cleaning also *~, *.x and gmon.out'
	@rm -f $(OBJECTS) *~ *.x gmon.out 
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

ode_io_utils.o: ode_io_utils.hh ode.hh model_ode.hh
ode_solve.o: ode_io_utils.hh ode.hh model_ode.hh 
ode.o: ode.hh 
model_ode.o: model_ode.hh ode.hh pa_macros.hh
