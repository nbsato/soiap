.SUFFIXES:
.SUFFIXES: .F .f90 .f .o

include make.inc

OBJS_o = algebra.o atomrelax_o.o input.o force_field.o keyvalue.o parms_o.o show_time.o structure_opt.o rfc5.o cif.o opt2_o_addtime.o

.f.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.f
.f90.o:
	$(FC) -c -o $@ $(FFLAGS_free) $*.f90
.F.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.F

all: opt2_o

opt2_o: $(OBJS_o)
	$(FC) -o $@  $^ $(LIBS)

atomrelax_o.o: parms_o.o 

force_field.o: parms_o.o algebra.o

input.o: keyvalue.o parms_o.o cif.o

structure_opt.o: parms_o.o 

clean:
	rm -f *.o *.mod 

cleanprog:
	rm -f opt2_o
