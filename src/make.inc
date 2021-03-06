FC = gfortran -fopenmp
FFLAGS = -O2 -cpp -g -fbacktrace
FFLAGS_fixed = $(FFLAGS) -ffixed-line-length-132
FFLAGS_free = $(FFLAGS) -ffree-line-length-256
LIBS = -llapack

# FC = ifort -qopenmp
# FFLAGS = -O2 -cpp -assume realloc_lhs -g -traceback
# FFLAGS_fixed = $(FFLAGS) -132
# FFLAGS_free = $(FFLAGS)
# LIBS = -mkl

# set if you use 'force_field 5'
# JMATGEN_ROOT = /path/to/jmatgen/root/directory
