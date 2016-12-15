import numpy as n
from numpy import (sqrt, linspace, asarray,
                    empty, hstack , vstack, pi)
from sys import argv

'''
Rectangular:
X, [:,0] = thru thickness
Y, [:,1] = Axial coord
Z, [:,2] = Angle coordinate theta (radians)
'''

ni3d_fine = n.load('./ConstructionFiles/ni3d_fine.npy')
ni3d_med = n.load('./ConstructionFiles/ni3d_med.npy')
ni3d_cors = n.load('./ConstructionFiles/ni3d_cors.npy')
ni3d_ref1_mid = n.load('./ConstructionFiles/ni3d_ref1_mid.npy')
ni3d_ref1_th = n.load('./ConstructionFiles/ni3d_ref1_th.npy')
ni3d_ref1_z = n.load('./ConstructionFiles/ni3d_ref1_z.npy')
ni3d_ref2_z = n.load('./ConstructionFiles/ni3d_ref2_z.npy')
nc_cors = n.load('./ConstructionFiles/node_coords_cors.npy')
nc_med = n.load('./ConstructionFiles/node_coords_med.npy')
nc_fine = n.load('./ConstructionFiles/node_coords_fine.npy')

try:
    worthless, Lg, Ltop, ODtop, ID, tg, R, num_el_fine_th, dt, eccen = argv
    Lg = float(Lg)  # Length of gage section (full, not considering symmetry plane)
    Ltop = float(Ltop)  # Length of thick section above radius/chamf
    ODtop = float(ODtop)    # Outer RADIUS of thick section
    ID = float(ID)  # Inner RADIUS
    tg = float(tg)  # thickness of test sxn
    R = float(R)    # Radius of chamf
    num_el_fine_th = int(num_el_fine_th)  # Num elements thru the test section thickness
except:
    pass

del nc_fine
del nc_med
del nc_cors

# Mesh definitions
num_el_fine_th = ni3d_fine.shape[0] - 1
num_el_fine_y = ni3d_fine.shape[1] - 1
num_el_fine_z = ni3d_fine.shape[2]
num_el_fine_tot = num_el_fine_z*num_el_fine_y*num_el_fine_th
num_el_med_th = ni3d_med.shape[0] - 1
num_el_med_y = ni3d_med.shape[1] - 1
num_el_med_z = ni3d_med.shape[2]
num_el_med_tot = num_el_med_z*num_el_med_y*num_el_med_th
num_el_cors_th = ni3d_cors.shape[0] - 1
num_el_cors_y = ni3d_cors.shape[1] - 1
num_el_cors_z = ni3d_cors.shape[2]
num_el_cors_tot = num_el_cors_z*num_el_cors_y*num_el_cors_th

########################################################    
################## Connnect on Fine ##################
########################################################
# Notice, the inner-most loop is "for i in num_el_fine_th",
# meaning I will connect thru-thickness first, then in z, then in y
elcon_fine = n.empty((num_el_fine_tot,8))
row = 0
for j in range(num_el_fine_y):
    for k in range(num_el_fine_z):
        for i in range(num_el_fine_th):
            index = n.array([[i,j+1,k], 
                        [i+1,j+1,k],
                        [i+1,j+1,k+1],
                        [i,j+1,k+1],
                        [i,j,k],
                        [i+1,j,k],
                        [i+1,j,k+1],
                        [i,j,k+1]]).astype(int)
            if k == num_el_fine_z -1:
                index[index[:,2] == k+1,2] = 0
            elcon_temp = [0,0,0,0,0,0,0,0]
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_fine[v[0],v[1],v[2]]
            elcon_fine[row] = elcon_temp
            row+=1
elnums = (n.arange(len(elcon_fine))+1)[:,None]
elcon_fine = hstack((elcon_fine,elnums)) 
print('Connected on Fine') 
########################################################    
################## Connnect on Med ##################
########################################################
elcon_med = n.empty((num_el_med_tot,8))
row = 0
for j in range(num_el_med_y):
    for k in range(num_el_med_z):
        for i in range(num_el_med_th):
            index = n.array([[i,j+1,k], 
                        [i+1,j+1,k],
                        [i+1,j+1,k+1],
                        [i,j+1,k+1],
                        [i,j,k],
                        [i+1,j,k],
                        [i+1,j,k+1],
                        [i,j,k+1]]).astype(int)
            elcon_temp = [0,0,0,0,0,0,0,0]
            if k == num_el_med_z -1:
                index[index[:,2] == k+1,2] = 0            
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_med[v[0],v[1],v[2]]
            elcon_med[row] = elcon_temp
            row+=1
elnums = (n.arange(len(elcon_med)) + n.max(elnums))[:,None]+1
elcon_med = hstack((elcon_med,elnums)) 
print('Connected on Med')
########################################################    
################## Connnect on Coarse ##################
########################################################
elcon_cors = n.empty((num_el_cors_tot,8))
row = 0
for j in (n.arange(num_el_cors_y)):
    for k in range(num_el_cors_z):
        for i in range(num_el_cors_th):
            index = n.array([[i,j+1,k], 
                        [i+1,j+1,k],
                        [i+1,j+1,k+1],
                        [i,j+1,k+1],
                        [i,j,k],
                        [i+1,j,k],
                        [i+1,j,k+1],
                        [i,j,k+1]]).astype(int)
            if k == num_el_cors_z -1:
                index[index[:,2] == k+1,2] = 0                        
            elcon_temp = [0,0,0,0,0,0,0,0]
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_cors[v[0],v[1],v[2]]
            elcon_cors[row] = elcon_temp
            row+=1
elnums = (n.arange(len(elcon_cors)) + n.max(elnums))[:,None]+1
elcon_cors = hstack((elcon_cors,elnums)) 
print('Connected on Coarse')
########################################################    
########## Connnect Ref1.  Fine to ref1_z ##############
########## and ref1_z to ref1_mid ######################
########################################################
# I have for i in... inside the if k%3 check so that I connect in
# complete thru-thickness stacks of each of the four similar-shape
# elements, then proceed circumferentially
j_fine = ni3d_fine.shape[1] - 1
j_mid = 0
j_ref = 0
elcon_ref1_z = n.empty((0,8))   # Too difficult to figure out the shape I need before hand 
for k in range(num_el_fine_z):
    if k%3 == 0:
        for i in range(num_el_fine_th): 
            nodes = [ni3d_ref1_mid[i,j_mid,k//3],
                     ni3d_ref1_mid[i+1,j_mid,k//3],
                     ni3d_ref1_z[i+1,j_ref,k+1],
                     ni3d_ref1_z[i,j_ref,k+1],
                     ni3d_fine[i,j_fine,k],
                     ni3d_fine[i+1,j_fine,k],
                     ni3d_fine[i+1,j_fine,k+1],
                     ni3d_fine[i,j_fine,k+1]
                    ]
            elcon_ref1_z = vstack((elcon_ref1_z,nodes))
    elif k%3 == 1:
        for i in range(num_el_fine_th):
            nodes = [ni3d_ref1_z[i,j_ref,k],
                     ni3d_ref1_z[i+1,j_ref,k],
                     ni3d_ref1_z[i+1,j_ref,k+1],
                     ni3d_ref1_z[i,j_ref,k+1],
                     ni3d_fine[i,j_fine,k],
                     ni3d_fine[i+1,j_fine,k],
                     ni3d_fine[i+1,j_fine,k+1],
                     ni3d_fine[i,j_fine,k+1]
                    ]
            elcon_ref1_z = vstack((elcon_ref1_z,nodes))
    elif k%3 == 2:
        for i in range(num_el_fine_th):
            # Exception for the last circumferential set
            if k == num_el_fine_z - 1:
                nodes = [ni3d_ref1_z[i,j_ref,k],
                         ni3d_ref1_z[i+1,j_ref,k],
                         ni3d_ref1_mid[i+1,j_mid,0],
                         ni3d_ref1_mid[i,j_mid,0],
                         ni3d_fine[i,j_fine,k],
                         ni3d_fine[i+1,j_fine,k],
                         ni3d_fine[i+1,j_fine,0],
                         ni3d_fine[i,j_fine,0]
                        ]
        
            else:   
                nodes = [ni3d_ref1_z[i,j_ref,k],
                         ni3d_ref1_z[i+1,j_ref,k],
                         ni3d_ref1_mid[i+1,j_mid,(k+1)//3],
                         ni3d_ref1_mid[i,j_mid,(k+1)//3],
                         ni3d_fine[i,j_fine,k],
                         ni3d_fine[i+1,j_fine,k],
                         ni3d_fine[i+1,j_fine,k+1],
                         ni3d_fine[i,j_fine,k+1]
                        ]
            elcon_ref1_z = vstack((elcon_ref1_z,nodes))
    # In a separate loop, connect ref1_z up to ref1_mid
    if k%3 == 0:
        for i in range(num_el_fine_th):            
            # Exception for the last the final z elements which need to 
            # connect back to the 0th znodes
            if k == num_el_fine_z - 3:
                nodes = [ni3d_ref1_mid[i,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,0],
                         ni3d_ref1_mid[i,j_mid,0],
                         ni3d_ref1_z[i,j_ref,k+1],
                         ni3d_ref1_z[i+1,j_ref,k+1],
                         ni3d_ref1_z[i+1,j_ref,k+2],
                         ni3d_ref1_z[i,j_ref,k+2]
                        ]
                elcon_ref1_z = vstack((elcon_ref1_z,nodes))
            else:
                nodes = [ni3d_ref1_mid[i,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,(k+3)//3],
                         ni3d_ref1_mid[i,j_mid,(k+3)//3],
                         ni3d_ref1_z[i,j_ref,k+1],
                         ni3d_ref1_z[i+1,j_ref,k+1],
                         ni3d_ref1_z[i+1,j_ref,k+2],
                         ni3d_ref1_z[i,j_ref,k+2]
                        ]
                elcon_ref1_z = vstack((elcon_ref1_z,nodes))    

elnums = (n.arange(len(elcon_ref1_z)) + n.max(elnums))[:,None]+1
elcon_ref1_z = hstack((elcon_ref1_z,elnums)) 

########################################################    
########## Connnect Ref1_mid to ref1_th ################
########## and ref1_th to med ######################
########################################################    
# Again, connect thru-thickness then proceed circumferentially
# Four loops for the four differently-shaped elements here
# 
j_med = 0
j_mid = 0
j_ref = 0
elcon_ref1_th = n.empty((0,8))
for k in range(num_el_med_z):
    for i in range(num_el_fine_th):
        # Exception for the last the final z elements which need to 
        # connect back to the 0th znodes
        if k == num_el_med_z - 1:
            k = -1 # We can do this since all elements have same z boundaries
        if i%3 == 0:
            nodes = [ni3d_med[i//3,j_med,k],
                     ni3d_ref1_th[i+1,j_ref,k],
                     ni3d_ref1_th[i+1,j_ref,k+1],
                     ni3d_med[i//3,j_med,k+1],
                     ni3d_ref1_mid[i,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k+1],
                     ni3d_ref1_mid[i,j_mid,k+1]
                    ]
            elcon_ref1_th = vstack((elcon_ref1_th,nodes))
        elif i%3 == 1:
            nodes = [ni3d_ref1_th[i,j_ref,k],
                     ni3d_ref1_th[i+1,j_ref,k],
                     ni3d_ref1_th[i+1,j_ref,k+1],
                     ni3d_ref1_th[i,j_ref,k+1],
                     ni3d_ref1_mid[i,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k+1],
                     ni3d_ref1_mid[i,j_mid,k+1]
                    ]
            elcon_ref1_th = vstack((elcon_ref1_th,nodes))
        elif i%3 == 2:
            nodes = [ni3d_ref1_th[i,j_ref,k],
                     ni3d_med[(i+1)//3,j_med,k],
                     ni3d_med[(i+1)//3,j_med,k+1],
                     ni3d_ref1_th[i,j_ref,k+1],
                     ni3d_ref1_mid[i,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k+1],
                     ni3d_ref1_mid[i,j_mid,k+1]
                    ]
            elcon_ref1_th = vstack((elcon_ref1_th,nodes))                
        if i%3 == 0:
            nodes = [ni3d_med[i//3,j_med,k],
                     ni3d_med[(i+3)//3,j_med,k],
                     ni3d_med[(i+3)//3,j_med,k+1],
                     ni3d_med[i//3,j_med,k+1],
                     ni3d_ref1_th[i+1,j_ref,k],
                     ni3d_ref1_th[i+2,j_ref,k],
                     ni3d_ref1_th[i+2,j_ref,k+1],
                     ni3d_ref1_th[i+1,j_ref,k+1]
                    ]
            elcon_ref1_th = vstack((elcon_ref1_th,nodes))
            
elnums = (n.arange(len(elcon_ref1_th)) + n.max(elnums))[:,None]+1
elcon_ref1_th = hstack((elcon_ref1_th,elnums))  
print('Connected on Ref1')
########################################################    
########## Connect med to ref2, ref2 to cors
########################################################
j_med = ni3d_med.shape[1] - 1
j_cors = 0
j_ref = 0
elcon_ref2 = n.empty((0,8))
for k in range(num_el_med_z):
    if k%3 == 0:
        for i in range(num_el_med_th): 
            nodes = [ni3d_cors[i,j_cors,k//3],
                     ni3d_cors[i+1,j_cors,k//3],
                     ni3d_ref2_z[i+1,j_ref,k+1],
                     ni3d_ref2_z[i,j_ref,k+1],
                     ni3d_med[i,j_med,k],
                     ni3d_med[i+1,j_med,k],
                     ni3d_med[i+1,j_med,k+1],
                     ni3d_med[i,j_med,k+1]
                    ]
            elcon_ref2 = vstack((elcon_ref2,nodes))
    elif k%3 == 1:
        for i in range(num_el_med_th):
            nodes = [ni3d_ref2_z[i,j_ref,k],
                     ni3d_ref2_z[i+1,j_ref,k],
                     ni3d_ref2_z[i+1,j_ref,k+1],
                     ni3d_ref2_z[i,j_ref,k+1],
                     ni3d_med[i,j_med,k],
                     ni3d_med[i+1,j_med,k],
                     ni3d_med[i+1,j_med,k+1],
                     ni3d_med[i,j_med,k+1]
                    ]
            elcon_ref2 = vstack((elcon_ref2,nodes))
    elif k%3 == 2:
        for i in range(num_el_med_th):
            # Exception for the last circumferential set
            if k == num_el_med_z - 1:
                nodes = [ni3d_ref2_z[i,j_ref,k],
                         ni3d_ref2_z[i+1,j_ref,k],
                         ni3d_cors[i+1,j_cors,0],
                         ni3d_cors[i,j_cors,0],
                         ni3d_med[i,j_med,k],
                         ni3d_med[i+1,j_med,k],
                         ni3d_med[i+1,j_med,0],
                         ni3d_med[i,j_med,0]
                        ]
        
            else:   
                nodes = [ni3d_ref2_z[i,j_ref,k],
                         ni3d_ref2_z[i+1,j_ref,k],
                         ni3d_cors[i+1,j_cors,(k+1)//3],
                         ni3d_cors[i,j_cors,(k+1)//3],
                         ni3d_med[i,j_med,k],
                         ni3d_med[i+1,j_med,k],
                         ni3d_med[i+1,j_med,k+1],
                         ni3d_med[i,j_med,k+1]
                        ]
            elcon_ref2 = vstack((elcon_ref2,nodes))
    # In a separate loop, connect ref2_z up to ref1_mid
    if k%3 == 0:
        for i in range(num_el_med_th):            
            if k == num_el_med_z - 3:
                nodes = [ni3d_cors[i,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,0],
                         ni3d_cors[i,j_cors,0],
                         ni3d_ref2_z[i,j_ref,k+1],
                         ni3d_ref2_z[i+1,j_ref,k+1],
                         ni3d_ref2_z[i+1,j_ref,k+2],
                         ni3d_ref2_z[i,j_ref,k+2]
                        ]
                elcon_ref2 = vstack((elcon_ref2,nodes))
            else:
                nodes = [ni3d_cors[i,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,(k+3)//3],
                         ni3d_cors[i,j_cors,(k+3)//3],
                         ni3d_ref2_z[i,j_ref,k+1],
                         ni3d_ref2_z[i+1,j_ref,k+1],
                         ni3d_ref2_z[i+1,j_ref,k+2],
                         ni3d_ref2_z[i,j_ref,k+2]
                        ]
                elcon_ref2 = vstack((elcon_ref2,nodes))    
elnums = (n.arange(len(elcon_ref2)) + n.max(elnums))[:,None]+1
elcon_ref2 = hstack((elcon_ref2,elnums)) 
print('Connected on Ref2')

elcon = n.vstack((elcon_fine,elcon_ref1_z,elcon_ref1_th,elcon_med,elcon_ref2,elcon_cors)).astype(int)
# save it, put the el numbers out front
# Reorder so that connectivitr processes in a different direction to get good mesh
elcon = elcon[:,[-1,3,2,1,0,7,6,5,4]]
#n.savetxt('Elements.dat',X = n.hstack((elcon[:,[-1]],elcon[:,:-1])) ,fmt='%.0f',delimiter=',')
#n.savetxt('elements_for_abaqus.dat',X = elcon ,fmt='%.0f',delimiter=',')
n.save('./ConstructionFiles/abaqus_elements.npy',elcon)


'''
def plotall():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    NC = n.load('./ConstructionFiles/node_coords_all.npy')
    X = NC[:,0]*n.cos(NC[:,2])
    Y = NC[:,0]*n.sin(NC[:,2])
    rgn = (NC[:,2]>=angle-1*dq_med) | (NC[:,2]<=1*dq_med)
    nodes = vstack((X[rgn],Y[rgn],NC[rgn,1],NC[rgn,-1])).T
    ax.plot(nodes[:,0],nodes[:,1],nodes[:,2],'b.',alpha=0.5)
    N = elcon_ref1_th[n.arange(-6,6),:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nodes[:,-1],i), nodes, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,1],NC[loc,2],'ko',alpha=0.3)
            line, = ax.plot(NC[loc,0],NC[loc,1],NC[loc,2],'ro',alpha=1)
            p.pause(0.5)
            line.remove()
    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.set_zlabel('Y')

def fineconcheck():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    node_coords_fine = n.load('./ConstructionFiles/node_coords_fine.npy')            
    rgn = (node_coords_fine[:,1] <=4*elht_fine) & (node_coords_fine[:,2] <=8*dq_fine)
    ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,2],node_coords_fine[rgn,1],'b.',alpha=0.4)
    N = elcon_fine[:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_fine[:,-1],i), node_coords_fine, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.3)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.25)
            line.remove()
    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.set_zlabel('Y')

def medcheckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    node_coords_med = n.load('./ConstructionFiles/node_coords_med.npy')            
    rgn = (node_coords_med[:,1] >=2*elht_med) & (node_coords_med[:,2] >=angle - 3*dq_med)
    ax.plot(node_coords_med[rgn,0],node_coords_med[rgn,2],node_coords_med[rgn,1],'b.',alpha=0.5)
    N = elcon_med[-6:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_med[:,-1],i), node_coords_med, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.3)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.25)
            line.remove()
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)

def corscheckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    ax.axis('off')
    node_coords_cors = n.load('./ConstructionFiles/node_coords_cors.npy')            
    rgn = (node_coords_cors[:,1] >=5*elht_cors) & (node_coords_cors[:,2] >=angle - 3*dq_cors)
    rgn = node_coords_cors[:,2] < angle
    ax.plot(node_coords_cors[rgn,0],node_coords_cors[rgn,2],node_coords_cors[rgn,1],'b.',alpha=0.3)
    N = elcon_cors[:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_cors[:,-1],i), node_coords_cors, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.05)
            line.remove()
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)

def ref1zcehckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    node_coords_all = n.load('./ConstructionFiles/node_coords_all.npy')    
    node_coords_fine = n.load('./ConstructionFiles/node_coords_fine.npy')            
    node_coords_ref1_th = n.load('./ConstructionFiles/node_coords_ref1_th.npy')            
    node_coords_ref1_mid = n.load('./ConstructionFiles/node_coords_ref1_mid.npy')            
    node_coords_ref1_z = n.load('./ConstructionFiles/node_coords_ref1_z.npy')      
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    rgn = ((node_indices_fine[:,1] == n.max(node_indices_fine[:,1])) &
            (node_indices_fine[:,2] <= 9) &
            (node_indices_fine[:,0] <= 6))
    zmax = n.max(node_coords_fine[rgn,2])
    xmax = n.max(node_coords_fine[rgn,0])
    ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,2],node_coords_fine[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_z[:,2] <= zmax)&
            (node_coords_ref1_z[:,0] <= xmax))
    ax.plot(node_coords_ref1_z[rgn,0],node_coords_ref1_z[rgn,2],node_coords_ref1_z[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_mid[:,2] <= zmax)&
            (node_coords_ref1_mid[:,0] <= xmax))
    ax.plot(node_coords_ref1_mid[rgn,0],node_coords_ref1_mid[rgn,2],node_coords_ref1_mid[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_th[:,2] <= zmax)&
            (node_coords_ref1_th[:,0] <= xmax))
    ax.plot(node_coords_ref1_th[rgn,0],node_coords_ref1_th[rgn,2],node_coords_ref1_th[rgn,1],'b.',alpha=0.5)
    N = elcon_ref1_z[4:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_all[:,-1],i), node_coords_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.5)
            line.remove()
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    
def ref1zcehckcon_end():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    node_coords_all = n.load('./ConstructionFiles/node_coords_all.npy')    
    node_coords_fine = n.load('./ConstructionFiles/node_coords_fine.npy')            
    node_coords_ref1_th = n.load('./ConstructionFiles/node_coords_ref1_th.npy')            
    node_coords_ref1_mid = n.load('./ConstructionFiles/node_coords_ref1_mid.npy')            
    node_coords_ref1_z = n.load('./ConstructionFiles/node_coords_ref1_z.npy')      
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    rgn = ((node_indices_fine[:,1] == n.max(node_indices_fine[:,1])))
    zmax = n.max(node_coords_fine[rgn,2])
    xmax = n.max(node_coords_fine[rgn,0])
    ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,2],node_coords_fine[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_z[:,2] <= zmax)&
            (node_coords_ref1_z[:,0] <= xmax))
    ax.plot(node_coords_ref1_z[rgn,0],node_coords_ref1_z[rgn,2],node_coords_ref1_z[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_mid[:,2] <= zmax)&
            (node_coords_ref1_mid[:,0] <= xmax))
    ax.plot(node_coords_ref1_mid[rgn,0],node_coords_ref1_mid[rgn,2],node_coords_ref1_mid[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_th[:,2] <= zmax)&
            (node_coords_ref1_th[:,0] <= xmax))
    ax.plot(node_coords_ref1_th[rgn,0],node_coords_ref1_th[rgn,2],node_coords_ref1_th[rgn,1],'b.',alpha=0.5)
    N = elcon_ref1_z[-4*num_el_fine_th:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_all[:,-1],i), node_coords_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.1)
            line.remove()
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)

def ref1thcheckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    node_coords_all = n.load('./ConstructionFiles/node_coords_all.npy')    
    node_coords_fine = n.load('./ConstructionFiles/node_coords_fine.npy')            
    node_coords_ref1_th = n.load('./ConstructionFiles/node_coords_ref1_th.npy')            
    node_coords_ref1_mid = n.load('./ConstructionFiles/node_coords_ref1_mid.npy')            
    node_coords_ref1_z = n.load('./ConstructionFiles/node_coords_ref1_z.npy')     
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    rgn = ((node_indices_fine[:,1] == n.max(node_indices_fine[:,1])) &
            (node_indices_fine[:,2] <= 9) &
            (node_indices_fine[:,0] <= 6))
    zmax = n.max(node_coords_fine[rgn,2])
    xmax = n.max(node_coords_fine[rgn,0])
    #ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,2],node_coords_fine[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_z[:,2] <= zmax)&
            (node_coords_ref1_z[:,0] <= xmax))
    ax.plot(node_coords_ref1_z[rgn,0],node_coords_ref1_z[rgn,2],node_coords_ref1_z[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_mid[:,2] <= zmax)&
            (node_coords_ref1_mid[:,0] <= xmax))
    ax.plot(node_coords_ref1_mid[rgn,0],node_coords_ref1_mid[rgn,2],node_coords_ref1_mid[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref1_th[:,2] <= zmax)&
            (node_coords_ref1_th[:,0] <= xmax))
    ax.plot(node_coords_ref1_th[rgn,0],node_coords_ref1_th[rgn,2],node_coords_ref1_th[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_med[:,2] <= zmax)&
            (node_coords_med[:,0] <= xmax)&
            (node_coords_med[:,1] == n.min(node_coords_med[:,1])))
    ax.plot(node_coords_med[rgn,0],node_coords_med[rgn,2],node_coords_med[rgn,1],'b.',alpha=0.5)    
    N = elcon_ref1_z[-8:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_all[:,-1],i), node_coords_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.5)
            line.remove()
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    
def ref2checkcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    node_coords_all = n.load('./ConstructionFiles/node_coords_all.npy')    
    node_indices_med = n.load('./ConstructionFiles/node_indices_med.npy')
    node_coords_cors = n.load('./ConstructionFiles/node_coords_cors.npy')                
    node_coords_ref2_z = n.load('./ConstructionFiles/node_coords_ref2_z.npy')                
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    rgn = ((node_indices_med[:,1] == n.max(node_indices_med[:,1])))
    zmax = n.max(node_coords_med[rgn,2])
    xmax = n.max(node_coords_med[rgn,0])
    ax.plot(node_coords_med[rgn,0],node_coords_med[rgn,2],node_coords_med[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_ref2_z[:,2] <= zmax)&
            (node_coords_ref2_z[:,0] <= xmax))
    ax.plot(node_coords_ref2_z[rgn,0],node_coords_ref2_z[rgn,2],node_coords_ref2_z[rgn,1],'b.',alpha=0.5)
    rgn = ((node_coords_cors[:,2] <= zmax)&
            (node_coords_cors[:,0] <= xmax)&
            (node_indices_cors[:,1] == 0))
    ax.plot(node_coords_cors[rgn,0],node_coords_cors[rgn,2],node_coords_cors[rgn,1],'b.',alpha=0.5)
    N = elcon_ref2[-4*num_el_med_th:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(node_coords_all[:,-1],i), node_coords_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.5)
            line.remove()
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)         
'''