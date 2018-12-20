.SUFFIXES:
.SUFFIXES: .F .f90 .f .o .a

include make.inc

ifeq ($(JMATGEN_ROOT),)
    LIBS_jmatgen = libjmatgen_dummy.a
else
    LIBS_jmatgen = $(JMATGEN_ROOT)/src/libjmatgen.a
    FFLAGS_jmatgen = -I$(JMATGEN_ROOT)/src
endif

OBJS_o = algebra.o atomrelax_o.o input.o force_field.o keyvalue.o parms_o.o show_time.o structure_opt.o rfc5.o cif.o soiap_o_addtime.o adp_base_module.o adp_kwu14_module.o periodic_lattice_module.o ewald.o lj_module.o

.f90.o:
	$(FC) -c -o $@ $(FFLAGS_free) $(FFLAGS_jmatgen) $*.f90
.F.o:
	$(FC) -c -o $@ $(FFLAGS_fixed) $*.F

all: soiap

soiap: $(OBJS_o) $(LIBS_jmatgen)
	$(FC) -o $@ $^ $(LIBS)

atomrelax_o.o: parms_o.o cif.o

force_field.o: parms_o.o algebra.o adp_kwu14_module.o periodic_lattice_module.o lj_module.o ewald.o $(LIBS_jmatgen)

input.o: keyvalue.o parms_o.o cif.o lj_module.o $(LIBS_jmatgen)

structure_opt.o: parms_o.o

parms_o.o: lj_module.o

adp_base_module.o: periodic_lattice_module.o

adp_kwu14_module.o: adp_base_module.o

lj_module.o: periodic_lattice_module.o

libjmatgen_dummy.a: jmatgen_dummy.o
	ar -r -v $@ $^

clean:
	rm -f *.o *.mod *.a

cleanprog:
	rm -f soiap
