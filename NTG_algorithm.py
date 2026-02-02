###################   Network-based tilling generation (NTG)   #########################
########################################################################################
#####
##### Generates a .fe file to be read by the surface evolver software consisting
##### in a 2D arrangement of cells with disorder in the location of the center of  
##### mass of each cell. Control on the density of cells p_holes and on the existence
##### of cell contacts in adjacent cells p_contact. Geometric disorder in the location 
##### of cells
#####
########################################################################################

########################################################################################
#####
##### Function: create_tiling_2D_Hexagon(L, p_contact, p_hole, name_output)
#####
##### L: length of the side of the square embedding the triangular/hexagonal cell 
#####   arangement in cell diameter units
##### p_contact: probability that two adjacent cells display an actual contact
##### p_hole: probability that a site is occupied by a cell
##### name_output: name of the generated file to be read by the surface evolver software
#####   Attention!! the extension must be ".fe"
#####
##### Output:  ".fe" file to be read by the evolver with a cell tiling with the
#####   desired properties 
#####
########################################################################################

import random
import numpy as np

def distance_2D(x,y):
    
    return np.sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]))

def Regular_Lattice_centers(L):

    #list of nodes

    list_Node_Geo=[]

    for i in range (0,L):
        if (i%2)==0:
            for k in range (0,L):
                list_Node_Geo.append([2*k, np.sqrt(3)*i])
        else:
            for k in range(1,L):
                list_Node_Geo.append([2*k-1, np.sqrt(3)*i])
    
    return list_Node_Geo

def Hexagons(centers, p_contact, p_hole):
    
    Geo_nodes=[]
    Hexagons_index=[]
    edges=[]

    for i in range(0, len(centers)):
        
        Hex_index=[]
        
        for k in range(0,6): ##### creating the nodes of the hexagon
            
            node=[]
            node=[centers[i][0]+np.cos(2*k*np.pi/6), centers[i][1]+np.sin(2*k*np.pi/6)]
            node_idx=len(Geo_nodes)+1
            
            T=1  #### check if the point exists!!
            
            for u in range(0, len(Geo_nodes)):
                
               if distance_2D(Geo_nodes[u][1], node)<0.001:
                   
                   if p_contact > random.random() :
                   
                       node_idx=u+1
                       T=0
                       break

            #### END check if the point exists!!

            if T==1: ### if the point does not exist, add its coordinates in the list      
                      
                Geo_nodes.append([node_idx, node])
                Hex_index.append(node_idx) #### add the node in the definition of the hexagon
            
            if T==0: ### if the node exists, only add it in the definition of the hexagon
                
                Hex_index.append(node_idx) 
        
        Hex_edges=[]

        for k in range(0,6): ##### creating the edges of the hexagon

            if k<5:
                edges.append([Hex_index[k],Hex_index[k+1]])
                Hex_edges.append(len(edges))
            
            if k==5:
                edges.append([Hex_index[k], Hex_index[0]])
                Hex_edges.append(len(edges))
        
        Hexagons_index.append(Hex_edges)
        
    for i in range(0, len(Geo_nodes)): ##randuniform: Geometric stochasticity
        
        Geo_nodes[i][1][0]=Geo_nodes[i][1][0]+random.uniform(-0.4,0.4)         
        Geo_nodes[i][1][1]=Geo_nodes[i][1][1]+random.uniform(-0.4,0.4)
        
    ############ removing cells at random according to the fraction CF    
    Hexagons_index_del=[]  
    
    for i in range(0, len(Hexagons_index)):
        
        if p_hole > random.random():
            Hexagons_index_del.append(Hexagons_index[i])
    
    ############ removing the edges from the cells that disappeared
    edges_del=[]
    
    for i in range(0, len(edges)):
        for k in range(0, len(Hexagons_index_del)):
            for u in range(0, len(Hexagons_index_del[k])):
                
                if Hexagons_index_del[k][u]==i+1:
                    edges_del.append([i+1, [edges[i][0], edges[i][1]]])
                    break
   
    ############# Cleaning the nodes that dropped after the deletion of cells        
    aux_edges=[]
    
    for i in range(0, len(edges_del)):        

        aux_edges.append(edges_del[i][1][0])
        aux_edges.append(edges_del[i][1][1])

    aux_edges=np.unique(aux_edges)  
    
    Geo_nodes_del=[]

    for i in range(0, len(aux_edges)):
        for k in range(0, len(Geo_nodes)):
            
            if aux_edges[i]==Geo_nodes[k][0]:
                Geo_nodes_del.append(Geo_nodes[k])
                del Geo_nodes[k]
                break
    
    return(Geo_nodes_del, edges_del, Hexagons_index_del)
    
### Print to a .fe file such that the surface evolver works!!
       
def print_hexagonal_tiling(Hexag, name_output):
    
    nodes=[]
    nodes=Hexag[0]
    
    File=open(name_output, 'w')
    File.write('STRING\n')
    File.write('space_dimension 2\n')
    File.write('\n')
    File.write('vertices\n')
    File.write('\n')
    
    for i in range(0, len(nodes)):
        
        File.write(str(nodes[i][0])+'  '+str(nodes[i][1][0])+' '+str(nodes[i][1][1]))
        File.write('\n')
           
    File.write('\n')
    File.write('\n')
    
    edges=[]
    edges=Hexag[1]
    File.write('edges  /* given by endpoints */\n')
    
    for i in range(0, len(edges)):

        File.write(str(edges[i][0])+'  '+str(edges[i][1][0])+' '+str(edges[i][1][1]))
        File.write('\n')
    
    File.write('\n')
    File.write('\n')
    File.write('faces  /* given by oriented edge loop */\n')
    
    Faces=[]
    Faces=Hexag[2]
    
    for i in range(0, len(Faces)):
        
        w=Faces[i][0]
        v=Faces[i][1]
        t=Faces[i][2]
        u=Faces[i][3]
        x=Faces[i][4]
        y=Faces[i][5]
        
        File.write(str(i+1)+'  '+str(w)+' '+str(v)+' '+str(t)+' '+str(u)+' '+str(x)+' '+str(y))
        File.write('\n')
    
    min_size=np.pi-0.2
    max_size=np.pi+0.25
    File.write('bodies  /* defined by their oriented faces */')
    File.write('\n')
    
    for i in range(0, len(Faces)):
        
        File.write(str(i+1)+'  '+str(i+1)+' '+'volume '+str(random.uniform(min_size,max_size)))
        File.write('\n')
    
    File.write('\n')
    File.write('read')
    File.close()
        
######################### Function: create tiling

def create_tiling_2D_Hexagon(L, p_contact, p_hole, name_output):

    g=[]
    g=Regular_Lattice_centers(L)
    Hexag=[]
    Hexag=Hexagons(g, p_contact, p_hole)
    print_hexagonal_tiling(Hexag, name_output)

######################### END

