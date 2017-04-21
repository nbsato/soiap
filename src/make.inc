
FC=gfortran  -fopenmp 
FFLAGS= -O2 -cpp -DNOSTOP  -g -DUSE_GEN 
FFLAGS_fixed= $(FFLAGS) -ffixed-line-length-132
FFLAGS_free= $(FFLAGS) -ffree-line-length-256
LIBS= -llapack -lblas 


#FC=ifort -openmp
#FFLAGS= -O2 -cpp -DUSE_GEN -DNOSTOP -132  -g -traceback 
#FFLAGS_fixed= $(FFLAGS) -132
#FFLAGS_free= $(FFLAGS) 
#LIBS= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

# mac
#FC=ifort 
#FFLAGS= -O0 -cpp -DUSE_GEN -DNOSTOP -132 -check all -g -traceback 
#FFLAGS_fixed= $(FFLAGS) -132
#FFLAGS_free= $(FFLAGS) 
#LIBS= -mkl
