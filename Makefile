FC 	= ifort
FF77	= ifort
#FC = gfortran
#FF77 = gfortran


HDIR 	= /gpfs/home/dwiczer/
#HDIR   = /home/david/

SOLDIR	= $(HDIR)Computation/DFBOLS/
SRCDIR  = $(HDIR)workspace/MW/


FDBGFLAGS = -g -mkl -init=snan -init=array 

FIFLAGS = -mkl -qopenmp -parallel -O3 -xhost 

FGFLAGS = -fopenmp -ffree-line-length-none 

GLIBS = -lblas -llapack -lgomp -L$(HDIR)Resources/lib -lnlopt
ILIBS = -L$(HDIR)Resources/lib -lnlopt

F77FLAGS= -c

SOLOBJS = altmov.o rescue_h.o update.o prelim_h.o bobyqa_h.o bobyqb_h.o trsbox_h.o

ifeq ($(FC),ifort)
	FCFLAGS = $(FIFLAGS)
	LIBS    = $(ILIBS)
else
	FCFLAGS = $(FGFLAGS)
	LIBS    = $(GLIBS)
endif

V0main: V0main.f90 $(SOLOBJS) V0para.o
	$(FC) $(FCFLAGS) V0main.f90 V0para.o $(SOLOBJS) -o V0main.out $(LIBS)

V0para.o: V0para.f90
	$(FC) -c -ffree-line-length-none $(SRCDIR)V0para.f90

# All the solver objects:
altmov.o :  $(SOLDIR)altmov.f
	$(FF77) $(F77FLAGS) $(SOLDIR)altmov.f

rescue_h.o :  $(SOLDIR)rescue_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)rescue_h.f

update.o :  $(SOLDIR)update.f
	$(FF77) $(F77FLAGS) $(SOLDIR)update.f
	
prelim_h.o :  $(SOLDIR)prelim_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)prelim_h.f	
	
bobyqa_h.o :  $(SOLDIR)bobyqa_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)bobyqa_h.f

bobyqb_h.o :  $(SOLDIR)bobyqb_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)bobyqb_h.f

trsbox_h.o :  $(SOLDIR)trsbox_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)trsbox_h.f
	
