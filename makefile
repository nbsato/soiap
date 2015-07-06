.SUFFIXES:
.SUFFIXES: .F .f90 .f .o

include make.inc

OBJS = parms.o keyvalue.o atomrelax_qmas.o

.f.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.f  -Iselect
.f90.o:
	$(FC) -c -o $@ $(FFLAGS_free) $*.f90  -Iselect 
.F.o:
	$(FC) -c -o $@ $(FFLAGS_free) $*.F  -Iselect 

all: opt

opt: opt.o $(OBJS)
	$(FC) -o $@  $^ $(LIBS)

opt.o: parms.o keyvalue.o atomrelax_qmas.o

atomrelax_qmas.o: parms.o atomrelax_qmas.o

parms.o: parms.f90

keyvalue.o: keyvalue.F

clean:
	rm -f *.o *.mod 

cleanprog:
	rm -f opt
