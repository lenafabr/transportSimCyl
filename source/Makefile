OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o particleutil.o transportutil.o
FC = gfortran
FLAGS = -O3 
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fcheck=all -Wno-tabs
LIBS =  -lblas -llapack
OSGFLAGS = -static
OSGLIBS =  -lblas -llapack -lgfortran -lm

# define directories
PDIR = ./profile
ODIR = .

# executable name
EXE = ../transportsimMT.exe
OSGEXE=../transportsimMT_osg.exe
TESTEXE = ../testrun.exe
PROFEXE = ../transportsimMT_prof.exe

.SUFFIXES: .o .f90 

# general suffix rules

$(PDIR)/%.o : tools/%.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -pg -c -o $@ $<

$(ODIR)/%.o : tools/%.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c -o $@ $<

$(PDIR)/%.o : %.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -pg -c -o $@ $<

$(ODIR)/%.o : %.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c -o $@ $<

default: $(addprefix $(ODIR)/,$(OBJS) main.o) # recipe for ODIR/OBJS
	$(FC) $(FLAGS) $(DEBUGFLAGS) $^ -o $(EXE) $(LIBS)

profile: $(addprefix $(PDIR)/,$(OBJS) main.o) # recipe for PDIR/OBJS
	$(FC) $(FLAGS) $(DEBUGFLAGS) -pg $^ -o $(PROFEXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(FLAGS) $(OSGFLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(OSGLIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod $(ODIR)/*.o $(ODIR)/*.mod $(PDIR)/*.o $(PDIR)/*.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
$(ODIR)/readkey.o: $(addprefix $(ODIR)/,inputparams.o keys.o genutil.o)
$(ODIR)/genutil.o: $(addprefix $(ODIR)/,mt19937.o)
$(ODIR)/transportutil.o: $(addprefix $(ODIR)/,keys.o mt19937.o particleutil.o)
$(ODIR)/particleutil.o: $(addprefix $(ODIR)/,keys.o genutil.o mt19937.o inputparams.o)

$(PDIR)/readkey.o: $(addprefix $(PDIR)/,inputparams.o keys.o genutil.o)
$(PDIR)/genutil.o: $(addprefix $(PDIR)/,mt19937.o)
$(PDIR)/transportutil.o: $(addprefix $(PDIR)/,keys.o mt19937.o particleutil.o)
$(PDIR)/particleutil.o: $(addprefix $(PDIR)/,keys.o genutil.o mt19937.o inputparams.o)