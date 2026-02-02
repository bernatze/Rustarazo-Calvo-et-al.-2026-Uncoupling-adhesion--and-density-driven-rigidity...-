Output of the scripts to be read by the _Surface evolver_ software:

Basic ref: Brakke, Kenneth A (1992), "The Surface Evolver", _Experimental Mathematics_, 1 (2): 141–165, doi:10.1080/10586458.1992.10504253
_______________________________

The scripts presented here have been used for simulations of tissues in the publication:

Rustarazo-Calvo, L. et al "Uncoupling adhesion- and density-driven rigidity transitions triggers epithelial organization in embryonic tissues" (2026) BiorXiv preprint: https://www.biorxiv.org/content/10.1101/2025.03.18.644006v1
_______________________________

## **Header script: NTG_algorithm.py**
 
### Network-based tilling generation (NTG)  

 Generates a .fe file to be read by the surface evolver software consisting in a 2D arrangement of cells with disorder in the location of the center of mass of each cell. Control on the density of cells p_holes and on the existence of cell contacts in adjacent cells p_contact. Geometric disorder in the location of cells

 Function: create_tiling_2D_Hexagon(L, p_contact, p_hole, name_output)

 L: length of the side of the square embedding the triangular/hexagonal cell 
   arangement in cell diameter units
   
 p_contact: probability that two adjacent cells display an actual contact
 
 p_hole: probability that a site is occupied by a cell
 
 name_output: name of the generated file to be read by the surface evolver software

--Attention!! the extension must be ".fe"

 Output: 
 
 ".fe" file to be read by the evolver with a cell tiling with the
   desired properties 
_______________________________

## **Header script: RDS_algrithm_Poly.py**

### Random Disk Spreading (RDS) algorithm (Polydisperse) 

 Function: create_random_tiling_2D(L, R, rho, delta, name_output)

 Generates a 2D random polydisperse tiling and writes it to a plain-text
  ".fe" file to be used as input for the Surface Evolver software.

 The generated tiling consists of cells with 3 posibles radii (polydispersity),
 randomly arranged inside a square domain with a prescribed target density.
 The generation algorithm enforces the target density exactly while allowing
 controlled geometric disorder.
 
 Parameters:
 
   L           : Linear size of the domain, expressed in number of cells
                 (the system contains approximately L^2 cells).
   
   R           : Reference radius used to set the mean cell size.
                 Individual cell radii are sampled around this value.
   
   rho         : Target packing density of the tiling.
   
   delta       : Maximum allowed deviation (in number of cells) from the
                 expected cells corresponding to the target density..
   
   name_output : Full path and name of the output file.

--Attention!! the extension must be ".fe"

 Output:
 
   A ".fe" file encoding a random 2D polydisperse tiling generated via the Random Disk Spreading (RDS) algorithm, readable by the Surface Evolver software.
_______________________________

  ## **Header script: RDS_algrithm.py**

  ### Random Disk Spreading (RDS) algorithm  (Monodisperse)

   Function: create_random_tiling_2D(L, R, rho, name_output)

 Generates a 2D random monodisperse tiling and writes it to a plain-text ".fe" file to be used as input for the Surface Evolver software.

 The generated tiling consists of circular cells with identical radius, arranged randomly inside a square domain with a prescribed target density. The generation algorithm enforces exact monodispersity and controlled geometric disorder.

 Parameters:

  L           : Linear size of the domain, expressed in number of cells
                 (the system contains approximately L^2 cells).

  R           : Radius of each cell (monodisperse configuration).
  
  rho         : Target packing density of the tiling.

 name_output : Full path and name of the output file.

 --Attention!! the extension must be ".fe"

 Output:

A ".fe" file encoding a random 2D monodisperse tiling generated via the Random Disk Spreading (RDS) algorithm, readable by the Surface Evolver software.
_______________________________

## **Header script: Surface_Evolver_simulation_instructions.py**

 This script generates plain-text instruction files used to run batches of _Surface Evolver_ simulations on a computing cluster.

 The generated files encode parameter sweeps over:

- Surface resolution (via the t parameter)
   
- Alpha values

- Tiling configurations


Output: .txt files to be consumed by Surface Evolver via command-line execution. (No evolver simulations are executed directly by this script).

