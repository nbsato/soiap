.SUFFIXES:
.SUFFIXES: .F .f90 .f .o

include make.inc

OBJS_o = parms_o.o keyvalue.o atomrelax_qmas_o.o show_time.o 
OBJS = parms.o keyvalue.o atomrelax_qmas.o show_time.o 

.f.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.f  -Iselect
.f90.o:
	$(FC) -c -o $@ $(FFLAGS_free) $*.f90  -Iselect 
.F.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.F  -Iselect 

all: opt2_o

opt2: opt2.o $(OBJS)
	$(FC) -o $@  $^ $(LIBS)

opt2_o: opt2_o.o $(OBJS_o)
	$(FC) -o $@  $^ $(LIBS)

opt: opt.o $(OBJS)
	$(FC) -o $@  $^ $(LIBS)

opt2.o: parms.o keyvalue.o atomrelax_qmas.o show_time.o  opt2_addtime.f90
	$(FC) -c -o $@  $(FFLAGS_free)  opt2_addtime.f90 
opt2_o.o: parms_o.o keyvalue.o atomrelax_qmas_o.o show_time.o  opt2_o_addtime.f90
	$(FC) -c -o $@  $(FFLAGS_free)  opt2_o_addtime.f90 
opt2_addtime.f90: opt2.f90
	awk -f converttool/addtime.awk opt2.f90 > opt2_addtime.f90
opt2_o_addtime.f90: opt2_o.f90
	awk -f converttool/addtime.awk opt2_o.f90 > opt2_o_addtime.f90

opt.o: parms.o keyvalue.o atomrelax_qmas.o show_time.o 

atomrelax_qmas.o: parms.o 

parms_o.o: parms_o.f90

keyvalue.o: keyvalue.F

clean:
	rm -f *.o *.mod 

cleanprog:
	rm -f opt
