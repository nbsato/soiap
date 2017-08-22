Instructions
============

Requirement
-----------

- Fortran compiler supporting the Fortran 2003 standard

Installation
------------

You can install **opt_cl2** by

1. cloning the Git repository,
2. switching to the branch `develop`,
3. moving in the directory `src`,
4. editing `make.inc`, and
5. running `make`.

An executable file `opt2_o` is created in the directory `src`.

Usage
-----

The following command executes **opt_cl2**:

~~~
opt2_o <input>
~~~

where `<input>` is the input file described below. Output files, which are also 
described below, are created in the current directory. Progress of the command 
is written in the standard output.

Input file
----------

Words after `!` are recognized as comments in the input file. Four samples are 
shown here. They are contained in the directory `work/samples`.

### Sample 1: manually given initial structure

~~~
unit_length 1 ! unit of lattice vectors and absolute atomic coords.
              ! 1: Ang or 2: Bohr

unit_vec_column 1 ! lattice vectors as 1: column or 2: raw vectors

unit_vec ! lattice vectors
( 2.5d0, 0    , 2.5d0 )
( 2.5d0, 2.5d0, 0     )
( 0    , 2.5d0, 2.5d0 )

volume_scale 1.0 ! scale the volume while keeping the shape of the cell

number_atom 2 ! number of atoms
number_element 1 ! number of elements
atom_pos 1 ! atomic positions in 1: relative or 2: absolute coords.

atom_list ! list of atoms
14  0    0    0     0 ! atomic number, x, y, z, fix
14  0.25 0.25 0.25  1 ! atomic number, x, y, z, relax

md_mode_cell 2 ! cell-relaxation method; 0: FIRE, 2: quenched MD, or 3: RFC5
number_max_relax_cell 100 ! max. number of the cell relaxation
mass_cell 17202.4d-2 ! cell "mass" in (electron mass)/Bohr

md_mode 4 ! atom-relaxation method
          ! 0: FIRE, 3: simple relax, or 4: quenched MD
number_max_relax 100 ! max. number of the atom relaxation
mass_scale 0 ! 0: do nothing or 1: rescale the atomic masses to that of C
max_displacement 0.25 ! max. displacement of atoms in Bohr

unit_time 1 ! unit of the time step; 1: fs or 2: Hartree a.u.
time_step 1.0 ! time step

external_stress_v 0.0 0.0 0.0 ! external pressure in GPa

th_force 5d-5 ! convergence threshold for the force in Hartree a.u.
th_stress 5d-7 ! convergence threshold for the stress in Hartree a.u.

force_field 1 ! force field
              ! 1: Stillinger-Weber for Si, 3: ZRL for Si-O-N-H, or 
              ! 4: ADP for Nd-Fe-B
~~~

### Sample 2: initial structure given by a CIF file

~~~
crystal initial.cif ! CIF file for the initial structure
symmetry 0 ! 0: not symmetrize displacements of the atoms or 1: symmetrize

md_mode_cell 2 ! cell-relaxation method; 0: FIRE, 2: quenched MD, or 3: RFC5
number_max_relax_cell 100 ! max. number of the cell relaxation
mass_cell 17202.4d-2 ! cell "mass" in (electron mass)/Bohr

md_mode 4 ! atom-relaxation method
          ! 0: FIRE, 3: simple relax, or 4: quenched MD
number_max_relax 100 ! max. number of the atom relaxation
mass_scale 0 ! 0: do nothing or 1: rescale the atomic masses to that of C
max_displacement 0.25 ! max. displacement of atoms in Bohr

unit_time 1 ! unit of the time step; 1: fs or 2: Hartree a.u.
time_step 1.0 ! time step

external_stress_v 0.0 0.0 0.0 ! external pressure in GPa

th_force 5d-5 ! convergence threshold for the force in Hartree a.u.
th_stress 5d-7 ! convergence threshold for the stress in Hartree a.u.

force_field 1 ! force field
              ! 1: Stillinger-Weber for Si, 3: ZRL for Si-O-N-H, or 
              ! 4: ADP for Nd-Fe-B
~~~

### Sample 3: variable-cell relaxation by the RFC5 method

The original RFC5 method relaxes the cell and atoms simultaneously. Within the 
implementation in **opt_cl2**, the cell and atoms are relaxed simultaneously 
only at the first atom-relaxation iteration of each cell-relaxation loop. 
Therefore, if you want to performe a variable-cell relaxation as is the original
RFC5 method, `number_max_relax` should be set to one. The cell and atoms are 
then relaxed simultaneously in all the iterations. See FAQ for details.

~~~
crystal initial.cif ! CIF file for the initial structure
symmetry 1 ! 0: not symmetrize displacements of the atoms or 1: symmetrize

md_mode_cell 3 ! cell-relaxation method; 0: FIRE, 2: quenched MD, or 3: RFC5
number_max_relax_cell 100 ! max. number of the cell relaxation
number_max_relax 1 ! max. number of the atom relaxation
max_displacement 0.1 ! max. displacement of atoms in Bohr

external_stress_v 0.0 0.0 0.0 ! external pressure in GPa

th_force 5d-5 ! convergence threshold for the force in Hartree a.u.
th_stress 5d-7 ! convergence threshold for the stress in Hartree a.u.

force_field 1 ! force field
              ! 1: Stillinger-Weber for Si, 3: ZRL for Si-O-N-H, or 
              ! 4: ADP for Nd-Fe-B
~~~

### Sample 4: fixed-cell relaxation by the RFC5 method

Within the implementation of the RFC5 method in **opt_cl2**, the cell and atoms 
are relaxed simultaneously at the first atom-relaxation iteration of each cell-
relaxation loop other than the first cell-relaxation iteration. Therefore, if 
`number_max_relax_cell` is set to one, the cell is not relaxed. See FAQ for 
details.

~~~
crystal initial.cif ! CIF file for the initial structure
symmetry 0 ! 0: not symmetrize displacements of the atoms or 1: symmetrize

md_mode_cell 3 ! cell-relaxation method; 0: FIRE, 2: quenched MD, or 3: RFC5
number_max_relax_cell 1 ! max. number of the cell relaxation
number_max_relax 100 ! max. number of the atom relaxation
max_displacement 0.1 ! max. displacement of atoms in Bohr

external_stress_v 0.0 0.0 0.0 ! external pressure in GPa

th_force 5d-5 ! convergence threshold for the force in Hartree a.u.
th_stress 5d-7 ! convergence threshold for the stress in Hartree a.u.

force_field 1 ! force field
              ! 1: Stillinger-Weber for Si, 3: ZRL for Si-O-N-H, or 
              ! 4: ADP for Nd-Fe-B
~~~

Output files
------------

There are four output files:

`log.struc`
:   Lattice vectors in Bohr and fractional atomic coordinates.

`log.tote`
:   Total energies in Hartree and volumes in Ang^3.

`log.frc`
:   Forces acting on atoms in Hartree a.u.

`log.strs`
:   Stress tensors in Hartree a.u.

In these files, `loopc` and `loopa` indicate the iteration numbers of cell and 
atom relaxation, respectively.

You can check convergence of a calculation in the standard output. The line 
starting with `QMD%frc converged.` means that the atom-relaxation iteration is 
converged, and the line starting with `QMD%strs converged.` means that the cell-
relaxation iteration is converged. If both iterations are converged, these two 
lines are written consecutively, and they are followed by three lines starting 
with `*** QMD%loopc =`, `tote,fmax,tote/natom,omega,omega/natom=`, and 
`!!!! normally end !!!!`. Summarizing the above, when both the cell- and atom-
relaxation iterations are converged, the standard output ends with, e.g., 

~~~
...
 QMD%loopc,QMD%loopa=          13           1  out of          100         100
show_time: tote_frc_strs:     0.001 sec
QMD%frc converged. QMD%loopc, QMD%loopa, QMD%fmax   13    1  0.9882E-16
QMD%strs converged. QMD%loopc, QMD%smax   13  0.8169E-09
 *** QMD%loopc =          13
tote,fmax,tote/natom,omega,omega/natom=       -0.3187202039        0.0000000000       -0.1593601019       40.0467640044       20.0233820022
 !!!! normally end !!!!
~~~

> **Note:**  
> The last line always says `normally end` even if it is not converged.

FAQ
---

### Which relaxation method should I use?

The RFC5 method is fastest, although it sometimes fails to converge. The 
quenched MD method is slower but safer than the RFC5 method.

### How large is the value of `mass_cell`?

Typically, the order of 1% of $M / \Omega_0^{1/3}$ is recommended, where $M$ is the 
total atomic mass per cell and $\Omega_0$ is the initial volume per cell. If 
`mass_cell` is too small, the structure usually results in isolated atoms and 
clusters. If `mass_cell` is too large, it takes long time to be converged.

### How can I perform a non-relaxation calculation?

If both `number_max_relax_cell` and `number_max_relax` are set to one, the cell 
and atoms are not relaxed. Quantities only of the initial structure are 
outputted. See the next topic for details.

### When does **opt_cl2** relax the cell and atoms?

For the RFC5 method,

| `loopc` | `loopa` | relax cell | relax atom |
|--------:|--------:|:----------:|:----------:|
|       1 |       1 |            |            |
|       1 |       2 |            |    relax   |
|       1 |       3 |            |    relax   |
|     ... |     ... |     ...    |     ...    |
|       1 |       M |            |    relax   |
|       2 |       1 |    relax   |    relax   |
|       2 |       2 |            |    relax   |
|       2 |       3 |            |    relax   |
|     ... |     ... |     ...    |     ...    |
|       2 |       N |            |    relax   |
|       3 |       1 |    relax   |    relax   |
|       3 |       2 |            |    relax   |
|       3 |       3 |            |    relax   |
|     ... |     ... |     ...    |     ...    |

For other methods,

| `loopc` | `loopa` | relax cell | relax atom |
|--------:|--------:|:----------:|:----------:|
|       1 |       1 |            |            |
|       1 |       2 |            |    relax   |
|       1 |       3 |            |    relax   |
|     ... |     ... |     ...    |     ...    |
|       1 |       M |            |    relax   |
|       2 |       1 |    relax   |            |
|       2 |       2 |            |    relax   |
|       2 |       3 |            |    relax   |
|     ... |     ... |     ...    |     ...    |
|       2 |       N |            |    relax   |
|       3 |       1 |    relax   |            |
|       3 |       2 |            |    relax   |
|       3 |       3 |            |    relax   |
|     ... |     ... |     ...    |     ...    |
