FC = gfortran -fopenmp
FFLAGS = -O2 -cpp -g -fbacktrace
FFLAGS_fixed = $(FFLAGS) -ffixed-line-length-132
FFLAGS_free = $(FFLAGS) -ffree-line-length-256
LIBS = -llapack

# FC = ifort -qopenmp
# FFLAGS = -O2 -cpp -g -traceback
# FFLAGS_fixed = $(FFLAGS) -132
# FFLAGS_free = $(FFLAGS)
# LIBS = -mkl

# mac
# FC = ifort
# FFLAGS = -O0 -cpp -check all -g -traceback
# FFLAGS_fixed = $(FFLAGS) -132
# FFLAGS_free = $(FFLAGS)
# LIBS = -mkl
