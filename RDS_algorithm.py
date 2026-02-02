##########     Random Disk Spreading (RDS) algorithm  (Monodisperse)  ##################
########################################################################################
#####
##### Function: create_random_tiling_2D(L, R, rho, name_output)
#####
##### Generates a 2D random monodisperse tiling and writes it to a plain-text
##### ".fe" file to be used as input for the Surface Evolver software.
#####
##### The generated tiling consists of circular cells with identical radius,
##### arranged randomly inside a square domain with a prescribed target density.
##### The generation algorithm enforces exact monodispersity and controlled
##### geometric disorder.
#####
##### Parameters:
#####   L           : Linear size of the domain, expressed in number of cells
#####                 (the system contains approximately L^2 cells).
#####   R           : Radius of each cell (monodisperse configuration).
#####   rho         : Target packing density of the tiling.
#####   name_output : Full path and name of the output file.
#####                 IMPORTANT: the file extension must be ".fe".
#####
##### Output:
#####   A ".fe" file encoding a random 2D monodisperse tiling generated via the
#####   Random Disk Spreading (RDS) algorithm, readable by the Surface Evolver software.
#####
########################################################################################


import random
import numpy as np

# =============================================================================
########### MONODISPERSE EMULSION ######################
# =============================================================================

# =============================================================================
# ############ Basic functions
# =============================================================================

def dist(u,v):
    
    return np.sqrt((u[0]-v[0])*(u[0]-v[0]) + (u[1]-v[1])*(u[1]-v[1]))

def angle(u,v):
    
    w=[]
    w=[v[0]-u[0],v[1]-u[1]]
    c=w[0]/np.sqrt(w[0]*w[0]+w[1]*w[1])
    if w[1]>=0:
        return(np.arccos(c))
        
    if w[1]<0:
        return(2*np.pi-np.arccos(c))
        
def midpoint(u,v):
    
    return [(u[0]+v[0])/2,(u[1]+v[1])/2]
    
# =============================================================================
# ########## Random spread of centers within the 2D plane
# =============================================================================

def generate_centers_2D(L,R,rho):

    Size=np.sqrt(3)*L/2
    N=int(rho*L*(L+1)) ####total number of disks accepted
    n=0
    Centers=[]
    volumes=[]
    count=100000
    vol=np.pi*0.21
    # 0.23
    
    c=0
        
    while N>n and c<count: ##### filling consecutively
        
        v=[]
        x=random.uniform(0.5, Size-0.5)
        y=random.uniform(0.5, Size-0.5)    
        v=[x,y]
    
        T=1
        
        for i in range(0,len(Centers)): ### NEGLIGIBLE VARIABILITY BETWEEN CELL SIZES
            
            # Optimal one 1.65 modified for super dense
            
            if dist(Centers[i],v)<1.54*R: ###### minimum dist between centers
                T=-1
                break
            
        
        if T>0:

            Centers.append(v)            
            volumes.append(0.99*vol) ###very little cell
                
            n=n+1
        c=c+1
        
        if c%(count/100)==0:
            print(c/(count/100))
    
    print(n)
    print(N)

# =============================================================================
#     ########### in the case the desired density has not been achived
#     ########### we search for holes that have not been randomly identified
# =============================================================================

    if n<N:
        
        print('checking potential holes')
        res=200
        coordinates=[] #### vector of the points 
        for i in range(0, res):
            
            for k in range(0,res):
              
                coordinates.append([0.5+(i/res)*(Size-1), 0.5+(k/res)*(Size-1)])
        
        random.shuffle(coordinates)
        
        for i in range(0, len(coordinates)): 
            
            d=[]
            
            for j in range(0, len(Centers)):
                d.append(dist(Centers[j],coordinates[i]))

            if min(d)>1.6*R:
                Centers.append(coordinates[i])
                volumes.append(0.99*vol)

                n=n+1
                
            if n>N:
                break
            
        print(n)
    return Centers, volumes

# =============================================================================
# #### returns the polygons with the right orientation
# #### Generates contacts if the distance between centers is <2R+epsilon
# =============================================================================

def create_disks(Centers, R): 
    
    res=15
    epsilon=0.05
    Polygons=[]
    points=[]
    points_Polygon=[]
    for k in range(0, len(Centers)):
        
        Polygons.append([])
        points_Polygon.append([])
        
    if len(Centers)>1:
        for k in range(0, len(Centers)):
            for i in range(k+1, len(Centers)):            
                if dist(Centers[k],Centers[i])<2*R+epsilon:
                    
                    mid=midpoint(Centers[k],Centers[i])
                    points.append(mid)                
                    Polygons[k].append(mid)
                    Polygons[i].append(mid)
                    points_Polygon[k].append(len(points)-1)
                    points_Polygon[i].append(len(points)-1)
                    
    for k in range(0, len(Centers)):
        for i in range(0,res):
            
            point=[]
            point=[Centers[k][0]+R*np.cos(2*i*np.pi/res),Centers[k][1]+R*np.sin(2*i*np.pi/res)]
            points.append(point)
            Polygons[k].append([point[0],point[1]])
            points_Polygon[k].append(len(points)-1)

    return Polygons, points_Polygon, points

# =============================================================================
# ####### topologically consistent orientations accounting 
# ####### for the contact points
# =============================================================================

def orientation(Centers, Polygons, points_Polygon, points):
    
    edges=[]
    oriented_faces=[]
    e=0
    
    for i in range(0, len(Centers)):
        
        angles_points=[]
        oriented_angles=[]
        
        for k in range(0, len(Polygons[i])): ##### create an ordered list of nodes to define links
            
            a=angle(Centers[i], Polygons[i][k])###### arccos of the point with respect center
            angles_points.append([a,points_Polygon[i][k]])

        oriented_angles.append(sorted(angles_points, key=lambda x:x[0]))       
        oriented_angles_single=np.transpose(oriented_angles[0])[1]
        oriented_angles_single=oriented_angles_single.tolist()
        oriented_angles_single.append(oriented_angles_single[0])
        
        ##### define the links of the face
        edges_face=[]
        for k in range(0, len(oriented_angles_single)-1):
            
            edge=[]
            edge=[oriented_angles_single[k]+1, oriented_angles_single[k+1]+1]
            ##### the +1 appears because the surface evolver starts at 1
            edges.append(edge)
            e=e+1
            edges_face.append(e) ### defining the edges of the face
        
        oriented_faces.append(edges_face)
        
    return points, edges,oriented_faces

# =============================================================================
# ### Print to a .fe file such that the surface evolver works!!
# =============================================================================

def print_random_tiling(nodes, edges, oriented_faces, volumes, name_output):
        
    File=open(name_output, 'w')
    File.write('STRING\n')
    File.write('space_dimension 2\n')
    File.write('\n')
    File.write('vertices\n')
    File.write('\n')
    
    for i in range(0, len(nodes)):
        
        File.write(str(i+1)+'  '+str(nodes[i][0])+' '+str(nodes[i][1]))
        File.write('\n')
           
    File.write('\n')
    File.write('\n')
    
    File.write('edges  /* given by endpoints */\n')
    
    for i in range(0, len(edges)):

        File.write(str(i+1)+'  '+str(int(edges[i][0]))+' '+str(int(edges[i][1])))
        File.write('\n')
    
    File.write('\n')
    File.write('\n')
    File.write('faces  /* given by oriented edge loop */\n')
    
    
    for i in range(0, len(oriented_faces)):
        
        File.write(str(i+1)) 
        
        for k in range (0, len(oriented_faces[i])):
            File.write(' '+str(oriented_faces[i][k]))

        File.write('\n')
    File.write('\n')

    File.write('bodies  /* defined by their oriented faces */')
    File.write('\n')
    
    for i in range(0, len(oriented_faces)):
        
        File.write(str(i+1)+'  '+str(i+1)+' '+'volume '+str(volumes[i]))
        File.write('\n')
    
    File.write('\n')
    File.write('read')
    File.close()
        
# =============================================================================
# ######################### Global Function: create tiling
# =============================================================================

def create_random_tiling_2D(L, R,rho, name_output):

    Centers=[]
    Tiling=[]
    nodes=[]
    edges=[]
    oriented_faces=[]
    Centers, volumes=generate_centers_2D(L,R,rho)
    Tiling=create_disks(Centers, R)
    nodes, edges,oriented_faces=orientation(Centers, Tiling[0], Tiling[1], Tiling[2])
    print_random_tiling(nodes, edges, oriented_faces, volumes, name_output)

######################### END ##########################################
    
# L=20
# R=0.5
# rho=0.8

# name_output='/Applications/Evolver270-OSX/fe/Lattices_FE/Fe_for_analysis_imaging/Underconstrained/Lattice_Cells_random_Til.fe'

# create_random_tiling_2D(L, R,rho, name_output)
