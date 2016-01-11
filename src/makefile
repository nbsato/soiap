.SUFFIXES:
.SUFFIXES: .F .f90 .f .o

include make.inc

OBJS_o = algebra.o atomrelax_qmas_o.o input.o force_field.o keyvalue.o parms_o.o show_time.o structure_opt.o

.f.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.f  -Iselect
.f90.o:
	$(FC) -c -o $@ $(FFLAGS_free) $*.f90  -Iselect 
.F.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.F  -Iselect 

all: opt2_o

opt2_o: opt2_o.o $(OBJS_o)
	$(FC) -o $@  $^ $(LIBS)

opt2_o.o: $(OBJS_o)  opt2_o_addtime.f90
	$(FC) -c -o $@  $(FFLAGS_free)  opt2_o_addtime.f90 
opt2_o_addtime.f90: opt2_o.f90
	awk -f converttool/addtime.awk opt2_o.f90 > opt2_o_addtime.f90

algebra.o: algebra.f90

atomrelax_qmas_o.o: parms_o.o 

force_field.o: parms_o.o algebra.o

input.o: keyvalue.o parms_o.o

keyvalue.o: keyvalue.F

parms_o.o: parms_o.f90

structure_opt.o: parms_o.o 

clean:
	rm -f *.o *.mod 

cleanprog:
	rm -f opt2_o
