# ------------------------------------------------------------------------
# Copyright (C) 2017 Nobuya Sato, Hiori Kino, and Takashi Miyake
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ------------------------------------------------------------------------

""" Soiap potential."""

# [note] set environmental variable SOIAP_COMMAND for soiap command

import os
import numpy as np

import ase
from ase.calculators.calculator import Calculator, all_changes

#---- calculation class for soiap
class Soiap(Calculator):

    implemented_properties = ["energy", "forces", "stress"]

    default_parameters = dict(
        number_max_relax=10,
        md_mode=4,
        number_max_relax_cell=0,
        external_stress_v=[0.0, 0.0, 0.0],
        unit_time=1,
        time_step=1.0,
        mass_scale=0,
        mass_cell=5.e+4,
        md_mode_cell=2,
        th_force=5.e-5,
        th_stress=5.e-7,
        force_field=1
    )

    #---- calculate enregy, forces, and stress
    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self,atoms,properties,system_changes)

        self.write_input()

        self.run()

        self.read_results()
    # end def calculate

    #---- write input file for the solve
    def write_input(self, atoms=None, properties=None, system_changes=None):
        # make convinience local variables
        cell      = self.atoms.cell
        symbols   = self.atoms.get_chemical_symbols()
        numbers   = self.atoms.get_atomic_numbers()
        positions = self.atoms.get_scaled_positions()
        constraint_indices = []
        for c in self.atoms.constraints:
            constraint_indices.extend(c.index)
        # end for

        numbers_uniq = np.unique(numbers)
        numbers_index = {}
        for e in range(len(numbers_uniq)):
            numbers_index[numbers_uniq[e]] = e+1
        # end for

        # write data to the input file
        fd = open( "input", "w" )
        fd.write("number_atom    %d\n" % len(self.atoms) )
        fd.write("number_element %d\n" % len(numbers_uniq) )
        fd.write("\n")

        fd.write("unit_vec\n")
        fd.write("( %f, %f, %f )\n" % ( cell[0][0], cell[0][1], cell[0][2]) )
        fd.write("( %f, %f, %f )\n" % ( cell[1][0], cell[1][1], cell[1][2]) )
        fd.write("( %f, %f, %f )\n" % ( cell[2][0], cell[2][1], cell[2][2]) )
        fd.write("\n")

        fd.write("unit_vec_column 1 ! column\n")
        fd.write("unit_length     1 ! Angstrom\n")
        fd.write("atom_pos        1 ! relative\n")
        fd.write("\n")

        fd.write("atom_list\n")
        for a in range(len(self.atoms)):
            fd.write("%2d %f %f %f %d %d ! %s\n" %
              (numbers[a], positions[a][0], positions[a][1], positions[a][2],
               numbers_index[numbers[a]],
               a not in constraint_indices, symbols[a] ) )
        # end for
        fd.write("\n")

        keys = sorted(self.parameters)
        for key in keys:
            if isinstance(self.parameters[key],list) :
                fd.write(
                  "%s " % key
                + "".join(" %r" %
                    self.parameters[key][i]
                      for i in range(len(self.parameters[key])))
                + "\n")
            else:
                fd.write("%s %r\n" % (key,self.parameters[key]) )
            # end if
        # end for

        fd.close()
    # end def write_input

    #---- run the solver
    def run(self):
        # exec the solver
        if "SOIAP_COMMAND" in os.environ:
            soiap = os.environ["SOIAP_COMMAND"]
            exitcode = os.system("%s input > log.stdout" % soiap )
        else:
            raise RuntimeError("Please set SOIAP_COMMAND environment variable")
        # end if

        if exitcode != 0:
            raise RuntimeError("Soiap exited with exit code: %d.  " % exitcode)
        # end if
    # end def run


    #---- read struct, enregy, forces, and stress
    def read_results(self):
        self.read_struct()
        self.read_energy()
        self.read_forces()
        self.read_stress()
    # end def read_results


    #---- read cell and positions
    def read_struct(self):
        # init local variable
        cell = [[0,0,0],[0,0,0],[0,0,0]]
        positions = []

        # get the latest values
        fd = open("log.struc","r")
        line = fd.readline()
        while line:
            if line.find("*** unit vectors [Bohr]:") != -1:
                positions = []

                line = fd.readline()
                terms = line.split()
                cell[0][0] = float(terms[0])
                cell[1][0] = float(terms[1])
                cell[2][0] = float(terms[2])

                line = fd.readline()
                terms = line.split()
                cell[0][1] = float(terms[0])
                cell[1][1] = float(terms[1])
                cell[2][1] = float(terms[2])

                line = fd.readline()
                terms = line.split()
                cell[0][2] = float(terms[0])
                cell[1][2] = float(terms[1])
                cell[2][2] = float(terms[2])

                line = fd.readline() # skip messeges ...
                for a in range(len(self.atoms)):
                    line = fd.readline()
                    terms = line.split()
                    position = [ float(terms[0]), float(terms[1]), float(terms[2]) ]
                    positions.append(position)
                # end for
            # end if
            line = fd.readline()
        # end while

        fd.close()

        # conversion of unit from bohr to angstrom
        for i in range(3):
           for j in range(3):
              self.atoms.cell[i][j] = cell[i][j]*ase.units.Bohr
           # end for
        # end for

        self.atoms.set_scaled_positions(positions)
    # end read_struct

    #---- read energy
    def read_energy(self):
        # init local variable
        energy = 0.0

        # get the latest values
        fd = open("log.tote","r")
        line = fd.readline()
        while line:
            terms = line.split()
            if len(terms) == 7:
                energy = float(terms[2])
            # end if
            line = fd.readline()
        # end while

        fd.close()

        # no conversion of unit from eV to eV
        self.results["energy"] = energy
    # end def read_energy

    #---- read forces
    def read_forces(self):
        # init local variable
        forces = []

        # get the latest values
        fd = open("log.frc","r")
        line = fd.readline()
        while line:
            if line.find("*** atomic forces: QMD%loopc") != -1:
                forces = []
                for a in range(len(self.atoms)):
                    line = fd.readline()
                    terms = line.split()
                    force = [ float(terms[0]), float(terms[1]), float(terms[2]) ]
                    forces.append( np.array(force) )
                # end for
            # end if
            line = fd.readline()
        # end while

        fd.close()

        # conversion of unit from eV/bohr to eV/angstrom
        for a in range(len(forces)):
           for j in range(3):
              forces[a][j] = forces[a][j]*(ase.units.eV/ase.units.Bohr)
           # end for
        # end for

        self.results["forces"] = np.array(forces)
    # end def read_forces

    #---- read stress
    def read_stress(self):
        # init local variable
        stress = np.zeros(6)

        # get the latest values
        fd = open("log.strs","r")
        line = fd.readline()
        while line:
            if line.find("QMD%loopc, QMD%smax") != -1:
                line = fd.readline()
                terms = line.split()
                stress[0] = float(terms[0]) # xx
                stress[5] = float(terms[1]) # xy
                stress[4] = float(terms[2]) # xz

                line = fd.readline()
                terms = line.split()
                stress[1] = float(terms[1]) # yy
                stress[3] = float(terms[2]) # yz

                line = fd.readline()
                terms = line.split()
                stress[2] = float(terms[2]) # zz
            # end if
            line = fd.readline()
        # end while

        fd.close()

        # conversion of unit from eV/bohr^3 to eV/angstrom^3
        for i in range(6):
            stress[i] = stress[i]*((-1.0)*ase.units.eV/ase.units.Bohr**3)
        # end for

        self.results["stress"] = stress
    # end def read_stress

# end class Soiap
