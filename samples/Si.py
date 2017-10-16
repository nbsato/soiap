from ase import Atoms
from ase.build import bulk
from ase.calculators.soiap import Soiap
from ase.constraints import FixAtoms

atoms = bulk('Si', 'diamond', a=6.429, cubic=True)
#constraint = []
#constraint.append( FixAtoms(indices=[0,1]) )
#constraint.append( FixAtoms(indices=[3,4]) )
#atoms.set_constraint( constraint )

calc = Soiap(
  number_max_relax=50,
  md_mode=0,
  number_max_relax_cell=100,
  external_stress_v=[0.0, 0.0, 0.0],
  unit_time=1,
  time_step=2.0,
  mass_scale=0,
  mass_cell=1.0e+2,
  md_mode_cell=3,
  th_force=1.e-6,
  th_stress=5.e-7,
  force_field=1
)
atoms.set_calculator(calc)

#---- calculate and show energy
energy = atoms.get_potential_energy()
print("energy[eV]\n %+f\n" % energy )

#---- show forces
forces = atoms.get_forces()
print("atomic forces[eV/angstrom]")
for a in range(len(forces)):
    print(" %+f, %+f, %+f" % (forces[a][0],forces[a][1],forces[a][2]))
# end for
print

#---- show stress
stress = atoms.get_stress()
print("stress[eV/angstrom^3]")
print(" %+f, %+f, %+f" % ( stress[0], stress[5], stress[4]) )
print(" %+f, %+f, %+f" % ( stress[5], stress[1], stress[3]) )
print(" %+f, %+f, %+f" % ( stress[4], stress[3], stress[2]) )
print

atoms_opt = calc.get_atoms()
#---- show cell
cell = atoms_opt.cell
print("unit_vec[angstrom]")
print(" %+f, %+f, %+f" % ( cell[0][0], cell[0][1], cell[0][2]) )
print(" %+f, %+f, %+f" % ( cell[1][0], cell[1][1], cell[1][2]) )
print(" %+f, %+f, %+f" % ( cell[2][0], cell[2][1], cell[2][2]) )
print

#---- show positions
coordinates = atoms_opt.get_scaled_positions()
print("internal coordinates")
for a in range(len(coordinates)):
    print(" %f, %f, %f" % 
      (coordinates[a][0],coordinates[a][1],coordinates[a][2]))
# end for
print
