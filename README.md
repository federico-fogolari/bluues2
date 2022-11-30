# Bluues2.0

Bluues2.0 performs most common analyses of the electrostatic properties of
molecules. 

The output of the program entails, depending on user’s requirement:
1) the generalized Born radius of each atom;
2) the electrostatic potential at the surface of the molecule mapped to solvent accessible atoms;
3) the solvent accessible surface in a PDB formatted file;
4) the electrostatic potential in a volume surrounding the molecule;
5) the electrostatic free energy and different contributions to it;
6) the pH-dependent properties (total charge and pH-dependent free energy of folding in the pH range -2 to 18;
7) the pKa of all ionizable groups;

The input file is a pqr file (or a PDB file with charges and radii in the two separated fields after coordinates), which can be generated in different ways:
for proteins and nucleic acids the program pdb2pqr by Nathan A Baker and collaborators (https://www.poissonboltzmann.org/) ; for small molecules openbabel can convert mol2 or sdf files in pqr; for general molecules softwares like Chimera or VMD can read structure and topology files and generate pqr files.

Bluues2 first generates the solvent accessible surface and then computes generalized Born radii, via a surface integral and then it uses generalized Born radii (using a finite radius test particle) to perform electrostic analyses. 
The generalized Born radii corresponding to the solvent accessible surface are extrapolated to those corresponding to the solvent excluded surface.

If the molecular (solvent excluded) surface is needed instead, the program msms should be installed and put in the path of the system. If the option -m followed by a basename for msms auxiliary generated files, the surfaces will be computed using msms.

CONTACT:  

Federico Fogolari  
Dipartimento di Scienze Matematiche, Informatiche e Fisiche  
Universita' di Udine  
Via delle Scienze 206  
33100 Udine - Italy  
Tel ++39 0432 494320    
E-mail federico.fogolari@uniud.it  


REFERENCE:

Please cite:  
?. ???,  ?. ???,  ?. ???,  ?. ??? and  ?. ???,   
Molecular electrostatics with the Generalized Born model.
A tutorial through examples with Bluues 2.0.
Molecules, ??, ??, 2022 

The methodology has been described in:

Fogolari F, Corazza A, Yarra V, Jalaru A, Viglino P and Esposito G.
Bluues: a program for the analysis of the electrostatic properties of proteins based on generalized Born radii
BMC Bioinformatics, 13 (4), 1-16, 2012.

The web-server version of the program is described in:

Walsh I, Minervini G, Corazza A, Esposito G, Tosatto S.C.E. and Fogolari F.
Bluues Server: electrostatic properties of wild-type and mutated protein structures. Bioinformatics. 28(16):2189-90. (2012)

============================================================================

COMPILATION:

The program is compiled with: 

> make bluues2

(Alternatively:
> cc bluues2.c -o bluues2 -lm -g -O3 
)

INSTALLATION:

> sudo make install

(Alternatively:
Just copy bluues2 in the local programs directory /usr/local/bin which should
have rwxr-xr-x permissions, otherwise just change them with

> sudo chmod 755 /usr/local/bin
> sudo cp bluues2 /usr/local/bin/bluues2

Just to be sure it can be run by anyone

> sudo chmod 755 /usr/local/bin/bluues2

Check that /usr/local/bin is in the path by

> echo $PATH

/usr/local/bin/ should be listed among others.
)

RUNNING BLUUES2 

After bluues is installed

> bluues2 

with no argument bluues2 provides a summary on how to run the program
and a list of the available options explained below:

Usage:
./bluues2 filename.pqr basename [Options]
Options:
-v (verbose mode)
-m filename.vert (msms .vert file with 3 lines header)
-pa x (msms area per point on SES, 0.1 default)
-p x (use a probe radius for surfacing of x, 1.5 A default)
-tcr x (use a test charge radius of x, 0.7 A default)
-mr x (minimum radius for surfacing of x, 1.0 A default)
-s x (use a salt radius of x, 2.0 default)
-pd x (inner dielectric constant, 4.0 default)
-sd x (outer dielectric constant, 78.54 default)
-kp x (factor for Still formula 4 default)
-i x (ionic strength (M), 0.15 M default)
-t x (temperature in K, 298.15 default)
-c x (cutoff for GBR6 calculation, 20.0 A default)
-c2 x (cutoff short for GBR6 calculation, 8.0 A default)
-g x (surface tension coefficient, 0.12 J/(A^2 mol) default)
-pka (compute pkas)
-pkadef filename (file containing pKa definitions)
-pkd x (inner dielectric constant for pka calculations, 20.0 default)
-pss x (number of MC steps per site in pka calculation, 1000)
-dx (output potential grid in dx format)
-nx x (number of grid points in x, 97 default)
-ny x (number of grid points in y, 97 default)
-nz x (number of grid points in z, 97 default)
-mesh x (grid mesh size, 1.0 default)
-srf (output surface points, potential, atom n. and normal vectors)
-srfpot (output average atomic surface potential)

USAGE EXAMPLES (input files and commands in the test directory):

1) COMPUTING GENERALIZED BORN RADII AND ELECTROSTATIC AND SOLVATION ENERGIES

Provided that a pqr file is available (e.g. protein.pqr)
the program is run in its basic form by specifying the input file
and the basename for output files:

> ./bluues2 protein.pqr protein_out

The command produces three output files:
protein_out.gbr -- the file containing the generalized Born (GB) 
                   radii for all atoms
protein_out.pqg -- a pdb file with occupancy and temperature factors fields
                   substituted by the charge and the computed GB radius
protein_out.solv_nrg -- a file containing the computed energy and the solvation energies of each atom (including self and Coulombic solvation energy)

The free energy components are: 
- the Coulombic energy computed assuming the inner dielectric constant filling all the space
- the sum of the Born energies of each atom
- the solvation screening energy of all Coulombic pairs
- the hydrophobic solvation energy computed as a surface tension
  gamma x solvent accessible surface area.

The sum of the latter three terms gives the total solvation energy.

The GB radii can be visualized in color using vmd (gbr.vmd contains a script for vmd) 

> vmd -e gbr.vmd 

This will create a picture (gbr.tga) with the radius coded in color red (1.00 A) to blue (10.0 A).
To keep the program interactive just remove the last command 'quit' from gbr.vmd.
 
2) COMPUTING THE POTENTIAL AT THE SURFACE

The potential at the surface of the protein is computed with the option -srfpot:
> ./bluues2 protein.pqr protein_out -srfpot

The output file protein_out.srfatpot is a pdb file which contains 
the surface potential in the temperature factor field. 
The surface potential is expressed as kJ/(mol * q).

The potential is computed using a test charge whose radius may be specified with the option -tcr followed by the value of the radius. The default for this value is 0.7 A.

We have chosen as an example the Antemnapaedia protein bound to DNA. 
The protein has different electrostatic properties on the face binding
DNA and on the opposite face.
This is visualized by 
> vmd -e srfpot.vmd 
which will create two pictures of the two faces mentioned above (srfpot_1.tga and srfpot_2.tga) with the potential coded in color red (-8.00 kJ/(mol * q)) 
to blue (8.00 kJ/(mol * q)).
The DNA which has been prepared in a different file, is read and shown in 
licorice representation.
To keep the program interactive just remove the last command 'quit' from srfpot.vmd.
 
3) COMPUTING THE POTENTIAL AROUND A MOLECULE 

The potential around a molecule is computed by the option -dx:

> ./bluues2 protein.pqr protein_out -dx

This will compute a 97x97x97 grid with mesh 1 A around the molecule and 
map the potential using the test charge as said above. 
The grid mesh size can be changed with the option -mesh
and the number of grid points in x, y and z, can be changed
by the options -nx, -ny, -nz respectively, each followed by the number of points in that dimension.

The map can be read and visualized. We load first the protein and display it showing its surface. Then we load the DNA and show it in licorice representation.
Finally we load the computed electrostatic map, and display a slice of it around the molecule. The command is:

> vmd -e heatmap.vmd

The position of the slice (axis and height) may by chosen by the user. Note that the map is computed in the absence of DNA which is shown to rationalize electrostatic properties.

Another possible visualization shows the isosurface curves around the molecule.
This is done with the command:

> vmd -e isosurface.vmd

4a) COMPUTING PKA SHIFTS FOR PROTEINS 

Bluues2 has built-in values for protein pKa.
Several quantities related with pKa are computed using the flag -pka

> ./bluues2 protein.pqr protein_out -pka

The output files are the following:

protein_out.pka - for each titratable site it reports:
- the computed pKa
- the unperturbed (reference) pKa in the isolated residue
- the shifts due to the difference in 
  -- self solvation energy
  -- interaction with partial charges in the molecule
  -- interaction with other titratable sites net charges
- the generalized Born radius of the titratable atom.

protein_out.ddg - for each pH value it reports
- the difference in electrostatic free energy of the molecule with
  respect to isolated titratable sites
- the overall average charge of the molecule.

protein_out.titration - for each titratable site it reports the 
 charge as a function of pH

4b) COMPUTING PKA SHIFTS FOR GENERAL MOLECULES

Bluues2 has the capability of performing pKa analysis for any molecule
for which the unperturbed pKa is provided by the user.
An example file with the data for proteins and DNA is distributed with the 
program.
pkadefs.txt contains the following information:
- the name of the residue
- the name of the atom representing the titratable group
- the unperturbed pKa, i.e. the pKa of the isolated residue
- a flag stating whether in the input pqr file the residue is ionized (1)
  or not (0)
- the charge of the residue when ionized
- the charge that must be added to neutralize the residue as provided
  in the input.pqr file (the last field is minus the product of the previous two  fields...)

As an example we compute the shift in pKa for the phosphates of DNA.

> ./bluues2 dna.pqr dna_out -pka -pkadef pkadefs.txt 

The computed shifts are positive and in the range of 0.6 to 1.5 pKa units.
The contribution of self solvation are small, being the phosphate exposed to solvent (the GB radius is small). Both background and site-site interactions 
shift the pKa to higher values. 
The details of the calculation (with the scaling applied) are reported in the reference papers.

