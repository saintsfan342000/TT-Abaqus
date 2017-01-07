import numpy as n
from numpy import (sqrt, linspace, asarray,
                    empty, hstack , vstack,
                    pi, meshgrid)
from sys import argv
import matplotlib.pyplot as p
from mpl_toolkits import mplot3d

'''
Node coordinates for 3D TT model
Specify Lg, Ltop, ODtop, ID, tg, numleth either on command line or in script.
Works by laying out the nodes in "orthogonal cylindrical" coordinates, 
then at the end using x=rcos(q),y-rsin(q) to map to a cylinder in cartesian space

In addition to a unique node number, each node gets an i,j,k index based on its position
in r, q, z space (referred to herein as x,z,y, respectively).

x or th, [:,0] = thru thickness, or radial coordinate
y, [:,1] = Axial coord
z or q, [:,2] = Angle coordinate theta
'''
try:
    worthless, Lg, Ltop, ODtop, ID, tg, R, num_el_fine_th, dt = argv
    Lg = float(Lg)  # Half-length of gage section
    Ltop = float(Ltop)  # Length of thick section above radius/chamf
    ODtop = float(ODtop)/2    # Outer RADIUS of thick section
    ID = float(ID)/2  # Inner RADIUS
    tg = float(tg)  # mean thickness of test sxn.  Must mean = 0.5*(tmin + tmax)
    R = float(R)    # Radius of chamf
    num_el_fine_th = int(num_el_fine_th)  # Num elements thru the test section thickness
    dt = float(dt) # Specify dt/tg.  
except:
    Lg = 0.4 / 2
    Ltop = 0.5  # Length of thick section above radius/chamf
    ODtop = 1.9685/2    # Radius of thick section
    ID = 1.75/2 # Inner Radius
    tg = .038
    R = .125    # Radius of chamf
    num_el_fine_th = 3 # Num elements thru the test section thicknes
    dt = 0

angle = 2*pi    # 2pi for full tube
coord_end_chamf = sqrt( R**2 - (ODtop-(ID+tg+R))**2 ) + Lg  # y-coord where chamfer reaches ODtop
ytop = coord_end_chamf + Ltop
start_ref1 = 2*Lg/3  # y-coord of last node in fine-mesh elements, where the refine starts
start_ref2 = coord_end_chamf  # y-coord of last node in med-mesh elements, where the refine starts

def CalcOD(Y):
    X = empty(Y.shape)
    rgn1 = (0<=Y) & (Y<=Lg)
    rgn2 = (Lg<Y) & (Y<=coord_end_chamf)
    rgn3 = (coord_end_chamf<Y)
    X[rgn1] = ID + tg
    X[rgn2] = -sqrt( R**2 - (Y[rgn2]-Lg)**2 ) + (ID+tg+R)
    X[rgn3] = ODtop
    return X

# Fine mesh definitions
elht_fine = tg / num_el_fine_th # Elemt height equal to element width
num_el_fine_y = int(start_ref1 // elht_fine + 1)   # Number of fine-mesh elements in y-dir
num_node_fine_y = num_el_fine_y + 1    # Num nodes in y-dir up to START of refine
num_node_fine_th = num_el_fine_th + 1
dq_fine = elht_fine/ID  # Angular step of an element
num_el_fine_z = int(angle//dq_fine)
if num_el_fine_z%9 != 0: 
    num_el_fine_z-=(num_el_fine_z%9)
dq_fine = angle/(num_el_fine_z-1)   # Corrected angular step
# medium mesh definitions
num_el_med_th = int(num_el_fine_th/3)
num_node_med_th = num_el_med_th + 1
elht_med = 3*elht_fine
start_med_y = start_ref1 + 2*elht_med
num_el_med_y = int((start_ref2 - start_med_y)//elht_med)
num_node_med_y = num_el_med_y + 1
dq_med=3*dq_fine
# medium mesh definitions
num_el_cors_th = num_el_med_th # Only refining height wise and circumfwise, not thru-thickness
num_node_cors_th = num_el_cors_th + 1
elht_cors = 3*elht_med
start_cors_y = start_ref2 + 1*elht_cors   # Again, not a thickness
num_el_cors_y = int((ytop - start_cors_y) // elht_cors)   # Number of elements up height in coarse rgn
num_node_cors_y = num_el_cors_y + 1   # Num nodes in y-dir from start of coarse to top
dq_cors = dq_med*3

# Layout the fine-mesh nodes
yspace = linspace(0, start_ref1, num_node_fine_y)   # Linspace of y-coordinates of fine nodes, up to start_ref1
yindspace = n.arange(len(yspace))
zspace = n.linspace(0, angle-dq_fine, num_el_fine_z)
zindspace = n.arange(len(zspace))
xindspace = n.arange(num_node_fine_th)
OD_fine = CalcOD(yspace)
node_coords_fine = empty((0,3))
node_indices_fine = empty((0,3),dtype=int)
for k,y in enumerate(yspace):
    # Node coords
    xo = ID+(.5*(dt*tg)*(1+n.cos(-pi*y/start_ref1))) # Here's the imperfection! 
    x = linspace(xo, OD_fine[k], num_node_fine_th)
    node_coords_fine = vstack(( node_coords_fine, vstack(map(n.ravel,n.meshgrid(x,y,zspace))).T ))
    # Node indices
    nod = asarray([yindspace[k]])
    node_indices_fine = vstack(( node_indices_fine, vstack(map(n.ravel,n.meshgrid(xindspace,k,zindspace))).T ))

nodenums = (n.arange(len(node_coords_fine))+1)[:,None]
node_coords_fine = hstack((node_coords_fine,nodenums))
node_indices_fine = hstack((node_indices_fine,nodenums))

# Now create the 3d-indexed array to store the node numbers
ni3d_fine = n.empty((len(xindspace),len(yindspace),len(zindspace)))
for v in node_indices_fine:
    ni3d_fine[v[0],v[1],v[2]] = v[3]

# Medium nodes    
yspace = linspace(start_med_y, start_ref2, num_node_med_y)
yindspace = n.arange(len(yspace))
zspace = zspace[::3] # Keep every third!
zindspace = n.arange(len(zspace))
xindspace = n.arange(num_node_med_th)
OD_med = CalcOD(yspace)
node_coords_med = empty((0,3))
node_indices_med = empty((0,3),dtype=int)
for k,y in enumerate(yspace):
    # Node coords
    x = linspace(ID, OD_med[k], num_node_med_th)
    node_coords_med = vstack(( node_coords_med, vstack(map(n.ravel,n.meshgrid(x,y,zspace))).T ))
    # Node indices
    nod = asarray([yindspace[k]])
    node_indices_med = vstack(( node_indices_med, vstack(map(n.ravel,n.meshgrid(xindspace,k,zindspace))).T ))

nodenums = (n.arange(len(node_coords_med)) + n.max(nodenums))[:,None]+1
node_coords_med = hstack((node_coords_med,nodenums))
node_indices_med = hstack((node_indices_med,nodenums))

ni3d_med = n.empty((len(xindspace),len(yindspace),len(zindspace)))
for v in node_indices_med:
    ni3d_med[v[0],v[1],v[2]] = v[3]

# coarse nodes    
yspace = linspace(start_cors_y, ytop, num_node_cors_y)   # Linspace of y-coordinates of cors nodes, up to start_ref1
yindspace = n.arange(len(yspace))
zspace = zspace[::3] # Keep every third!
zindspace = n.arange(len(zspace))
xindspace = n.arange(num_node_cors_th)
OD_cors = CalcOD(yspace)
node_coords_cors = empty((0,3))
node_indices_cors = empty((0,3),dtype=int)
for k,y in enumerate(yspace):
    # Node coords
    x = linspace(ID, OD_cors[k], num_node_cors_th)
    node_coords_cors = vstack(( node_coords_cors, vstack(map(n.ravel,n.meshgrid(x,y,zspace))).T ))
    # Node indices
    nod = asarray([yindspace[k]])
    node_indices_cors = vstack(( node_indices_cors, vstack(map(n.ravel,n.meshgrid(xindspace,k,zindspace))).T ))

nodenums = (n.arange(len(node_coords_cors)) + n.max(nodenums))[:,None]+1
node_coords_cors = hstack((node_coords_cors,nodenums))
node_indices_cors = hstack((node_indices_cors,nodenums))        

ni3d_cors = n.empty((len(xindspace),len(yindspace),len(zindspace)))
for v in node_indices_cors:
    ni3d_cors[v[0],v[1],v[2]] = v[3]
    
'''
for k in 1:
    ## View all nodes at q = 0, looking down the circumferential die
    %matplotlib
    p.figure()
    rgn = NI[:,2]==0
    p.plot(NC[rgn,0],NC[rgn,1],'.')
    for k,i in enumerate(NC[rgn,0]):
        p.text(i,NC[rgn][k,1],'({:.0f},{:.0f})'.format(NI[rgn][k,0],NI[rgn][k,1]),ha='center',va='center',alpha=0.4)

    # Plot all thru-thickness nodes but just a small, fast-to-plot rgn of the other two directions
    # and write the node's indices 
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    rgn = (node_indices_fine[:,1] <= 1)&(node_indices_fine[:,2]<=3)
    ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,2],node_coords_fine[rgn,1],'b.',alpha=0.6)
    C = node_indices_fine[rgn].copy()
    for k,i in enumerate(node_coords_fine[rgn]):
        ax.text(i[0],i[2],i[1],'{},{},{}'.format(C[k,0],C[k,2],C[k,1]))
        
    # Demonstrate ref1
    p.close('all')
    %matplotlib
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    rgn = ((node_indices_fine[:,1] == n.max(node_indices_fine[:,1])) &
            (node_indices_fine[:,2] <= 9) &
            (node_indices_fine[:,0] <= 3))
    zmax = n.max(node_coords_fine[rgn,2])
    xmax = n.max(node_coords_fine[rgn,0])
    ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,2],node_coords_fine[rgn,1],'b.',alpha=0.5)
    for k,i in enumerate(node_coords_fine[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_fine[rgn][k,0],
                                            node_indices_fine[rgn][k,1],
                                            node_indices_fine[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='b',fontsize=10)
    rgn = ((node_coords_ref1_z[:,2] <= zmax)&
            (node_coords_ref1_z[:,0] <= xmax))
    ax.plot(node_coords_ref1_z[rgn,0],node_coords_ref1_z[rgn,2],node_coords_ref1_z[rgn,1],'r.',alpha=0.5)
    for k,i in enumerate(node_coords_ref1_z[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_ref1_z[rgn][k,0],
                                            node_indices_ref1_z[rgn][k,1],
                                            node_indices_ref1_z[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='r',fontsize=10)
    rgn = ((node_coords_ref1_mid[:,2] <= zmax)&
            (node_coords_ref1_mid[:,0] <= xmax))
    ax.plot(node_coords_ref1_mid[rgn,0],node_coords_ref1_mid[rgn,2],node_coords_ref1_mid[rgn,1],'g.',alpha=0.5)
    for k,i in enumerate(node_coords_ref1_mid[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_ref1_mid[rgn][k,0],
                                            node_indices_ref1_mid[rgn][k,1],
                                            node_indices_ref1_mid[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='g',fontsize=10)
    rgn = ((node_coords_ref1_th[:,2] <= zmax)&
            (node_coords_ref1_th[:,0] <= xmax))
    ax.plot(node_coords_ref1_th[rgn,0],node_coords_ref1_th[rgn,2],node_coords_ref1_th[rgn,1],'k.',alpha=0.5)
    for k,i in enumerate(node_coords_ref1_th[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_ref1_th[rgn][k,0],
                                            node_indices_ref1_th[rgn][k,1],
                                            node_indices_ref1_th[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='k',fontsize=10)
    rgn = ((node_coords_med[:,2] <= zmax)&
            (node_coords_med[:,0] <= xmax)&
            (node_indices_med[:,1] == n.min(node_indices_med[:,1])))
    ax.plot(node_coords_med[rgn,0],node_coords_med[rgn,2],node_coords_med[rgn,1],'.',alpha=1,color='orange')
    for k,i in enumerate(node_coords_med[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_med[rgn][k,0],
                                            node_indices_med[rgn][k,1],
                                            node_indices_med[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=1,color='orange',fontsize=10)

    # Demonstrate ref2
    p.close('all')
    %matplotlib
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Z',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    rgn = ((node_indices_med[:,1] == n.max(node_indices_med[:,1])) &
            (node_indices_med[:,2] <= 9) &
            (node_indices_med[:,0] <= 3))
    zmax = n.max(node_coords_med[rgn,2])
    xmax = n.max(node_coords_med[rgn,0])
    ax.plot(node_coords_med[rgn,0],node_coords_med[rgn,2],node_coords_med[rgn,1],'b.',alpha=0.5)
    for k,i in enumerate(node_coords_med[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_med[rgn][k,0],
                                            node_indices_med[rgn][k,1],
                                            node_indices_med[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='b',fontsize=10)
    rgn = ((node_coords_ref2_z[:,2] <= zmax)&
            (node_coords_ref2_z[:,0] <= xmax))
    ax.plot(node_coords_ref2_z[rgn,0],node_coords_ref2_z[rgn,2],node_coords_ref2_z[rgn,1],'r.',alpha=0.5)
    for k,i in enumerate(node_coords_ref2_z[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_ref2_z[rgn][k,0],
                                            node_indices_ref2_z[rgn][k,1],
                                            node_indices_ref2_z[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='r',fontsize=10)
    rgn = ((node_coords_cors[:,2] <= zmax)&
            (node_coords_cors[:,0] <= xmax)&
            (node_indices_cors[:,1] == n.min(node_indices_cors[:,1])))
    ax.plot(node_coords_cors[rgn,0],node_coords_cors[rgn,2],node_coords_cors[rgn,1],'.',alpha=1,color='g')
    for k,i in enumerate(node_coords_cors[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(node_indices_cors[rgn][k,0],
                                            node_indices_cors[rgn][k,1],
                                            node_indices_cors[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=1,color='g',fontsize=10)        
        
        
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('Y',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)
    rgn = (node_indices_fine[:,1] <= 1)&(node_indices_fine[:,2]<=3)
    ax.plot(node_coords_fine[rgn,0],node_coords_fine[rgn,1],node_coords_fine[rgn,2],'b.',alpha=0.6)
    C = node_indices_fine[rgn].copy()
    for k,i in enumerate(node_coords_fine[rgn]):
        ax.text(i[0],i[1],i[2],'{},{},{}'.format(C[k,0],C[k,1],C[k,2]))
    
'''

# First refine (start_ref1) circumferential (z) nodes
y = start_ref1 + elht_med/2
# locs are those fine nodes that are at the top y-coord of the fine section
locs = node_coords_fine[:,1] == n.max(node_coords_fine[:,1])
x = n.sort( n.unique( node_coords_fine[locs,0] )) 
z = n.sort( n.unique( node_coords_fine[locs,2] ))
node_coords_ref1_z = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Indices
y = n.array([0])
x = n.arange(len(x))
z = n.arange(len(z))
node_indices_ref1_z = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Now take out node with yindex = 0,3,6,9,...
node_coords_ref1_z = node_coords_ref1_z[ node_indices_ref1_z[:,2]%3 != 0, :]
node_indices_ref1_z = node_indices_ref1_z[ node_indices_ref1_z[:,2]%3 != 0, :]
# nodenums
nodenums = (n.arange(len(node_coords_ref1_z)) + n.max(nodenums))[:,None]+1
node_coords_ref1_z = hstack((node_coords_ref1_z,nodenums))
node_indices_ref1_z = hstack((node_indices_ref1_z,nodenums))        
#3d index array
ni3d_ref1_z = n.empty((len(x),len(y),len(z)))
for v in node_indices_ref1_z:
    ni3d_ref1_z[v[0],v[1],v[2]] = v[3]

# First refine (start_ref1)  mid-level nodes (btwn the circumf. and the thickness refine)
# thickness (x) like fine, circumf. (z) like med
y = start_ref1 + elht_med
x = n.sort( n.unique( node_coords_fine[locs,0] ))  # locs defined above to account for imperfection
z = n.sort( n.unique( node_coords_med[:,2] ))
node_coords_ref1_mid = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Indices
y = n.array([0])
x = n.arange(len(x))
z = n.arange(len(z))
node_indices_ref1_mid = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
#nodenums
nodenums = (n.arange(len(node_coords_ref1_mid)) + n.max(nodenums))[:,None]+1
node_coords_ref1_mid = hstack((node_coords_ref1_mid,nodenums))
node_indices_ref1_mid = hstack((node_indices_ref1_mid,nodenums))
#3d index array
ni3d_ref1_mid = n.empty((len(x),len(y),len(z)))
for v in node_indices_ref1_mid:
    ni3d_ref1_mid[v[0],v[1],v[2]] = v[3]

# First refine (start_ref1) thickness (th) nodes
y = start_ref1 + 3*elht_med/2
x = n.sort( n.unique( node_coords_fine[locs,0] )) # locs defined above to account for imperfection
z = n.sort( n.unique( node_coords_med[node_indices_med[:,1]==0,2] ))
node_coords_ref1_th = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Indices
y = n.array([0])
x = n.arange(len(x))
z = n.arange(len(z))
node_indices_ref1_th = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Now take out node with yindex = 0,3,6,9,...
node_coords_ref1_th = node_coords_ref1_th[ node_indices_ref1_th[:,0]%3 != 0, :]
node_indices_ref1_th = node_indices_ref1_th[ node_indices_ref1_th[:,0]%3 != 0, :]
#nodenums
nodenums = (n.arange(len(node_coords_ref1_th)) + n.max(nodenums))[:,None]+1
node_coords_ref1_th = hstack((node_coords_ref1_th,nodenums))
node_indices_ref1_th = hstack((node_indices_ref1_th,nodenums))     
#3d index array
ni3d_ref1_th = n.empty((len(x),len(y),len(z)))
for v in node_indices_ref1_th:
    ni3d_ref1_th[v[0],v[1],v[2]] = v[3]


# Second refine (start_ref2) circumferential (z) nodes
y = start_ref2 + elht_cors/2
# med crosses the radius from test sxn to ends, so unique won't work unless we
# also specify that we want the nodes at the top of the med section
x = n.sort( n.unique( node_coords_med[node_coords_med[:,1]==start_ref2,0] ))
z = n.sort( n.unique( node_coords_med[node_coords_med[:,1]==start_ref2,2] ))
node_coords_ref2_z = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Indices
y = n.array([0])
x = n.arange(len(x))
z = n.arange(len(z))
node_indices_ref2_z = vstack(map(n.ravel,n.meshgrid(x,y,z))).T
# Now take out node with yindex = 0,3,6,9,...
node_coords_ref2_z = node_coords_ref2_z[ node_indices_ref2_z[:,2]%3 != 0, :]
node_indices_ref2_z = node_indices_ref2_z[ node_indices_ref2_z[:,2]%3 != 0, :]
#nodenums
nodenums = (n.arange(len(node_coords_ref2_z)) + n.max(nodenums))[:,None]+1
node_coords_ref2_z = hstack((node_coords_ref2_z,nodenums))
node_indices_ref2_z = hstack((node_indices_ref2_z,nodenums))   
#3d index array
ni3d_ref2_z = n.empty((len(x),len(y),len(z)))
for v in node_indices_ref2_z:
    ni3d_ref2_z[v[0],v[1],v[2]] = v[3]

# All node coords
NC = vstack((node_coords_fine,node_coords_med,node_coords_cors,node_coords_ref1_z,node_coords_ref1_mid,node_coords_ref1_th,node_coords_ref2_z))
NI = vstack((node_indices_fine,node_indices_med,node_indices_cors,node_indices_ref1_z,node_indices_ref1_mid,node_indices_ref1_th,node_indices_ref2_z))

'''
for k in [1]:
    # Plot from ID looking out (x == 0)
    rgn = NI[:,0] == 0
    p.plot(NC[rgn,2],NC[rgn,1],'.')
    rgn = (node_coords_ref1_z[:,0] == ID)
    p.plot(node_coords_ref1_z[rgn][:,2],node_coords_ref1_z[rgn][:,1],'r.')
    rgn = (node_coords_ref1_mid[:,0] == ID)
    p.plot(node_coords_ref1_mid[rgn][:,2],node_coords_ref1_mid[rgn][:,1],'g.')
    rgn = (node_indices_ref1_th[:,0] == 1)
    p.plot(node_coords_ref1_th[rgn][:,2],node_coords_ref1_th[rgn][:,1],'yo')
    rgn = (node_coords_ref2_z[:,0] == ID)
    p.plot(node_coords_ref2_z[rgn][:,2],node_coords_ref2_z[rgn][:,1],'r.')
    
    # Plot looking in hoop direction
    p.figure()
    rgn = NI[:,2] == 0
    p.plot(NC[rgn,0],NC[rgn,1],'.')
    rgn = (node_indices_ref1_z[:,2] == 1)
    p.plot(node_coords_ref1_z[rgn][:,0],node_coords_ref1_z[rgn][:,1],'r.')
    rgn = (node_indices_ref1_mid[:,2] == 0)
    p.plot(node_coords_ref1_mid[rgn][:,0],node_coords_ref1_mid[rgn][:,1],'g.')
    rgn = (node_indices_ref1_th[:,2] == 0)
    p.plot(node_coords_ref1_th[rgn][:,0],node_coords_ref1_th[rgn][:,1],'yo')
    rgn = (node_indices_ref2_z[:,2] == 1)
    p.plot(node_coords_ref2_z[rgn][:,0],node_coords_ref2_z[rgn][:,1],'r.')

    # 3D Polar
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    rgn = NC[:,2] <= pi/2
    c = NC[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'.',alpha=0.15)
    rgn = node_coords_ref1_z[:,2] <= pi/2
    c = node_coords_ref1_z[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'r.',alpha=0.5)
    rgn = node_coords_ref1_mid[:,2] <= pi/2
    c = node_coords_ref1_mid[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'g.',alpha=0.5)
    rgn = node_coords_ref1_th[:,2] <= pi/2
    c = node_coords_ref1_th[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'k.',alpha=0.5)
    rgn = node_coords_ref2_z[:,2] <= pi/2
    c = node_coords_ref2_z[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'r.',alpha=0.5)

   
'''

'''
NC_polar_X = NC[:,0]*n.cos(NC[:,2])
NC_polar_Y = NC[:,0]*n.sin(NC[:,2])
fig = p.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(NC_polar_X,NC_polar_Y,NC[:,1],'.',alpha=0.3)
ax.set_xlabel('x')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()

p.figure()
p.plot(NC_polar_X,NC_polar_Y,'.',alpha=0.3)
p.axis('equal')
'''

nodelists = ['node_coords_cors', 'node_coords_fine', 'node_coords_med',
            'node_coords_ref1_mid', 'node_coords_ref1_th', 'node_coords_ref1_z',
            'node_coords_ref2_z', 'node_indices_cors', 'node_indices_fine',
            'node_indices_med', 'node_indices_ref1_mid', 'node_indices_ref1_th',
            'node_indices_ref1_z', 'node_indices_ref2_z', 'NC','NI',
            'ni3d_fine', 'ni3d_med', 'ni3d_cors', 'ni3d_ref1_mid', 
            'ni3d_ref1_th', 'ni3d_ref1_z', 'ni3d_ref2_z','node_abaqus']
 
 
for k, F in enumerate(nodelists):
    if F.rfind('coords') == 5:
        #n.savetxt('./ConstructionFiles/{}.dat'.format(F),X=eval(F),fmt = '%.12f',delimiter=',')
        n.save('./ConstructionFiles/{}'.format(F),eval(F))
    elif F.rfind('NC') == 0:
        #n.savetxt('./ConstructionFiles/{}.dat'.format('node_coords_all'),X=eval(F),fmt = '%.12f,%.12f,%.12f,%.0f')
        n.save('./ConstructionFiles/{}'.format('node_coords_all'),eval(F)) 
    elif F.rfind('NI') == 0:
        #n.savetxt('./ConstructionFiles/{}.dat'.format('node_indices_all'),X=eval(F),fmt = '%.0f,%.0f,%.0f,%.0f')
        n.save('./ConstructionFiles/{}'.format('node_indices_all'),eval(F).astype(int)) 
    elif F.rfind('ni3d_') == 0:
        n.save('./ConstructionFiles/{}'.format(F),eval(F).astype(int)) 
    elif F.rfind('abaqus') == 5:
        '''
        # node list for abaqus input file
        X = NC[:,0]*n.cos(NC[:,2])  
        Y = NC[:,0]*n.sin(NC[:,2])
        data = n.vstack((NC[:,-1],X,Y,NC[:,1])).T
        n.savetxt('nodes.dat',X=data,fmt='%.0f, %.12f, %.12f, %.12f')        
        '''
        pass
    elif F.rfind('indices') == 5:
        #n.savetxt('./ConstructionFiles/{}.dat'.format(F),X=eval(F),fmt = '%.0f',delimiter=',')
        n.save('./ConstructionFiles/{}'.format(F),eval(F).astype(int)) 
        pass    
    else:
        raise('Something bad in nodelist')
