from ase.calculators.emt import EMT
from ase.build import bulk

atoms = bulk('C', 'diamond', a=6.429, cubic=True)

calc = EMT()
atoms.set_calculator(calc)

energy = atoms.get_potential_energy()
print("energy\n %+f\n" % energy )
