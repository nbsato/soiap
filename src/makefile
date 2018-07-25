.SUFFIXES:
.SUFFIXES: .F .f90 .f .o

include make.inc

PATH_jmatgen = ../../jmatgen/src
LIB_jmatgen =  -L$(PATH_jmatgen) -ljmatgen
FFLAGS_jmatgen = -I$(PATH_jmatgen)

OBJS_o = algebra.o atomrelax_o.o input.o force_field.o keyvalue.o parms_o.o show_time.o structure_opt.o rfc5.o cif.o soiap_o_addtime.o adp_base_module.o adp_kwu14_module.o periodic_lattice_module.o

.f90.o:
	$(FC) -c -o $@ $(FFLAGS_free) $(FFLAGS_jmatgen) $*.f90
.F.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.F

all: soiap

soiap: $(OBJS_o)
	$(FC) -o $@ $^ $(LIB_jmatgen) $(LIBS)

atomrelax_o.o: parms_o.o cif.o

force_field.o: parms_o.o algebra.o adp_kwu14_module.o periodic_lattice_module.o $(PATH_jmatgen)/libjmatgen.a

input.o: keyvalue.o parms_o.o cif.o

structure_opt.o: parms_o.o 

adp_base_module.o: periodic_lattice_module.o

adp_kwu14_module.o: adp_base_module.o

clean:
	rm -f *.o *.mod 

cleanprog:
	rm -f soiap
