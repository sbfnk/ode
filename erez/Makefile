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

SOURCES     = main.cpp ode.cpp io_utils.cpp van_der_pol.cpp
OBJECTS     = $(SOURCES:.cpp=.o)
EXEC        = out

#--------------------------------------------------------------------

all:  $(EXEC)

$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@.x
	@echo

%.o : %.cpp
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) *.cpp

clean:
	@echo ... cleaning $(OBJECTS)
	@rm -f $(OBJECTS)
	@echo ... done

cleaner:
	@echo '... cleaning also *~ and gmon.out'
	@rm -f $(OBJECTS) *~ gmon.out
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

io_utils.o: io_utils.h ode.h van_der_pol.h
main.o: io_utils.h ode.h van_der_pol.h  
ode.o: ode.h 
van_der_pol: van_der_pol.h ode.h
