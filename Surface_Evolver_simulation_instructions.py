####################   Surface_Evolver_simulation_instructions.py   ############################
########################################################################################
#####
##### This script generates plain-text instruction files used to run batches of
##### Surface Evolver simulations on a computing cluster.
#####
##### The generated files encode parameter sweeps over:
#####   - Surface resolution (via the t parameter)
#####   - Alpha values
#####   - Tiling configurations
#####
##### Output: .txt files to be consumed by Surface Evolver via command-line execution.
##### (No evolver simulations are executed directly by this script).
#####
########################################################################################


#### ------------------- Simulation parameters -------------------- ###

#### Input/output paths ####
PathIn  = 'folder_path_in'  # Directory containing the original .fe files and where Evolver instruction files will be save.
PathOut = 'folder_path_out'  # Directory where simulation instructions will be save.

#### Alpha values to be simulated ####
alphas =['0.700','0.725','0.750','0.775','0.800','0.825','0.850','0.860','0.864','0.868','0.872','0.880','0.890','0.905','0.920','0.935','0.950','0.975']

#### ---------------- Tiling parameters  ---------------- ###

N_tot = 50  # Total number of tilings to simulate

# Expected filename format:
# 'Lattice_L{L}Cells_R{R}_rho{rho}_N{i}.fe'

L = 26      # Size of the network (~ L^2 cells)
R = 0.5     # Disk radius (code assumes RDS generation at R = 0.5)
rho = [0.8] # Target tiling density (list structure to allow sweeps)


#############################################################################

#### ---------------- Evolver resolution parameters ---------------- ####

# The surface resolution in Evolver is controlled by the t parameter (Smaller t values correspond to higher surface resolution).

# In our simulation, we use two different approaches: fixed and variable t.

# 1) Fixed resolution (constant t value):
#    - Set t_fix = True
#    - Only t_min is used as a constant t value (desired resolution)

# 2) Variable resolution (t sampled from a Weibull distribution):
#    - Set t_fix = False
#    - t values are sampled between t_min and t_max
#    - Weibull parameters (shape, scale) must be provided

# --- Mode selection ---
t_fix = True  # True: fixed t value | False: Weibull-distributed t values

# --- Fixed t configuration ---
t_min = 0.005   # Minimum allowed t value (highest surface resolution)
# --- Variable t configuration (used only if t_fix = False) ---
t_max = 0.14    # Maximum allowed t value (lowest surface resolution)
# --- Weibull distribution parameters (used only if t_fix = False) ---
shape = 1.4     # Weibull shape parameter
scale = 0.035   # Weibull scale parameter

# Label used to track the selected t-resolution range in filenames and outputs
t_label = 'minp005maxp14'  # example corresponds to t_min = 0.005 and t_max = 0.14
it = 20  # Number of iterations per alpha value


############## ----------------- CODE ---------------------- ################

import random
import numpy as np

p_hole_string=[]
for m in rho:
    p_hole_string.append(str(m))

bins=np.linspace(t_min,t_max,n_bins)

#Simulation generator begins here:
N=-1
for gg in range (0,N_tot):
    N = N + 1  
    for k in range (0,len(p_hole_string)):    

        MasterFileName = 'Simulator_L{}_Frac{}_Alphas_N{}_t{}_It{}.txt'.format(L,p_hole_string[k] ,N,t_label,it)
        
        File=open(PathIn+MasterFileName, 'w', newline='\n')
        
        #Use the following three lines only to run the visual snapshots
        
        #File.write('echo "s\n') #instructions to be passed to Evolver as input begin here
        #File.write('q\n') #evolver command
        #File.write('o\n') #evolver command
        File.write("echo 'o\n") #evolver command
        File.write('r\n') #evolver command to break segments
        File.write('r\n')
        File.write('r\n')

        # To stabilize the simulation, we initialize Evolver at the highest segment
        # resolution used (t_min) and run four iterations at alpha = 1.0 (hard spheres) 
        # before exploring different alpha or t-resolution values.

        File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*1.0)\n')
        File.write(f't {t_min}; g 100; o; g 100; r; o; g 100\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*1.0)\n')
        File.write(f't {t_min}; g 100; o; g 100; r; o; g 100\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*1.0)\n')
        File.write(f't {t_min}; g 100; o; g 100; r; o; g 100\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
        File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*1.0)\n')
        File.write(f't {t_min}; g 100; o; g 100; r; o; g 100\n')

        # Export the files corresponding to the stabilized initial tiling.
        
        # This saves the .fe file after the stabilization conditions are applied.                 
        filename_fe ='Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}.fe'.format(L,R,p_hole_string[k], '1', N, 4,t_label)
        File.write('DUMP "{}{}"\n'.format(PathOut, filename_fe))
        
        # Export the simulation data to a .txt file.
        # (For storage reasons, this step can be skipped by commenting it out if not needed.)
        # filename = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}.txt'.format(L,R,p_hole_string[k], '1', N, 4,t_label)   
        # File.write('foreach edges ff where sum(ff.faces,1) == 2 do list ff.facets >> "{}{}"\n'.format(PathOut, filename))
        
        # Export the simulation data to a .ps file.
        # (For storage reasons, this step can be skipped by commenting it out if not needed.)                   
        # filename_ps = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}.ps'.format(L,R,p_hole_string[k], '1', N, 4,t_label)         
        # File.write('FULL_BOUNDING_BOX ON; ps_stringwidth:=0.00004; POSTSCRIPT "{}{}"\n'.format(PathOut, filename_ps))
        

        # For the TCJ analysis, these lines export data to make quantitative analysis easier. (For storage reasons, these 3 steps can be skipped by commenting it out if not needed.) 
        # This exports the network contacts.
        filename_txt = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}_Contacts.txt'.format(L,R,p_hole_string[k], '1', N, 4,t_label) 
        comand='{foreach ff.facet gg do print gg.id}'
        File.write('foreach edges ff where sum(ff.faces,1) == 2 do {} >> "{}{}"\n'.format(comand,PathOut, filename_txt))
        # This export the number of cells
        filename_Ncells = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}_Ncells.txt'.format(L,R,p_hole_string[k], '1', N, 4,t_label) 
        File.write('foreach facet ff do print ff.id >> "{}{}"\n'.format(PathOut, filename_Ncells))
        # This export the number of TCJ
        filename_TCJ = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}_TCJ.txt'.format(L,R,p_hole_string[k], '1', N, 4,t_label) 
        comand_TCJ='{foreach vv.facet gg do print gg.id}'
        File.write('foreach vertices vv where sum(vv.facets,1) == 3 do {} >> "{}{}"\n'.format(comand_TCJ,PathOut, filename_TCJ))
        

        # Start main simulations at varying alpha values after stabilization
        for aa in reversed(alphas):
                
            for l in range (0,it-1): # Run 'it-1' iterations at the desired alpha for different resolutions.

                if t_fix==True:

                    # Evolver commands corresponding to the chosen t resolution values.
                    File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
                    File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*{})\n'.format(aa))
                    File.write(f't {t_min}; g 100; o; g 100; r; o; g 100\n')        
                    l = l+1 
                    
                else:
                    # The Weibull distribution is used to generate t resolution values constrained between t_min and t_max.
                    n_rand=np.random.weibull(shape,1)*scale
                    bin_index=np.digitize(n_rand+t_min/2,bins)-1
                    t_rand=np.round(bins[bin_index],3)
    
                    # Evolver commands corresponding to the chosen t resolution values.
                    File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
                    File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*{})\n'.format(aa))
                    File.write(f't { t_rand[0]}; g 100; o; g 100; r; o; g 100\n')        
                    l = l+1 

            # Run one additional iteration at the highest resolution (t_min)
            # to ensure the final resolution is consistent across all alphas and tilings.
            File.write('foreach edges ff where sum(ff.faces,1) == 1 do ff.tension:=1.0\n')
            File.write('foreach edges ff where sum(ff.faces,1) == 2 do ff.tension:=(2*{})\n'.format(aa))
            File.write(f't { t_min}; g 100; o; g 100; r; o; g 100\n')

            # Export the files corresponding to the simulated alpha.
            
            # The resulting .fe file reflects the configuration after the simulation conditions.
            filename_fe ='Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}.fe'.format(L,R,p_hole_string[k], aa, N, l,t_label)
            File.write('DUMP "{}{}"\n'.format(PathOut, filename_fe))

            # Export the simulation data to a .txt file.
            # (For storage reasons, this step can be skipped by commenting it out if not needed.)
            # filename = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}.txt'.format(L,R,p_hole_string[k], aa, N, l,t_label)   
            # File.write('foreach edges ff where sum(ff.faces,1) == 2 do list ff.facets >> "{}{}"\n'.format(PathOut, filename))

            # Export the simulation data to a .ps file.
            # (For storage reasons, this step can be skipped by commenting it out if not needed.)     
            # filename_ps = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}.ps'.format(L,R,p_hole_string[k], aa, N, l,t_label)                    
            # File.write('FULL_BOUNDING_BOX ON; ps_stringwidth:=0.00004; POSTSCRIPT "{}{}"\n'.format(PathOut, filename_ps))

            # For the TCJ analysis, these lines export data to make quantitative analysis easier. (For storage reasons, these 3 steps can be skipped by commenting it out if not needed.)
            # This exports the network contacts.
            filename_txt = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}_Contacts.txt'.format(L,R,p_hole_string[k], aa, N, l,t_label) 
            comand='{foreach ff.facet gg do print gg.id}'
            File.write('foreach edges ff where sum(ff.faces,1) == 2 do {} >> "{}{}"\n'.format(comand,PathOut, filename_txt))
            # This export the number of cells
            filename_Ncells = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}_Ncells.txt'.format(L,R,p_hole_string[k], aa, N, l,t_label) 
            File.write('foreach facet ff do print ff.id >> "{}{}"\n'.format(PathOut, filename_Ncells))
            # This export the number of TCJ
            filename_TCJ = 'Fcells_L{}_R{}_rho{}_Alpha{}_N{}_Iteration{}_t{}_TCJ.txt'.format(L,R,p_hole_string[k], aa, N, l,t_label) 
            comand_TCJ='{foreach vv.facet gg do print gg.id}'
            File.write('foreach vertices vv where sum(vv.facets,1) == 3 do {} >> "{}{}"\n'.format(comand_TCJ,PathOut, filename_TCJ))
        

        File.write('q\n') # evolver command
        
        # Expected filename format:
        # 'Lattice_L{L}Cells_R{R}_rho{rho}_N{i}.fe'.format(L, R, rho[k], i)
        
        # Adjust the path to the directory where Evolver is installed on the cluster
        inputname = '/usr/local/bin/evolver {}Lattice_L{}Cells_R{}_rho{}_N{}.fe'.format(PathIn,L,R,p_hole_string[k],N)    
        
        File.write("q' | {}\n".format(inputname)) #evolver command
                      
    
        File.write('\n')
        # File.write('read')
        File.close()


############## ----------------- CODE ---------------------- ################

#############################################################################