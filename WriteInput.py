import numpy as n
from numpy import (linspace, asarray, in1d, empty,
                    hstack , vstack, pi, compress)
from sys import argv
from pandas import read_excel
import time

'''
Writes the input file
'''

if len(argv) == 8:
    expt = int(argv[1])
    inpname = argv[2]
    constit = argv[3]
    num_el_fine_th = argv[4]
    dt = argv[5]
    ecc = argv[6]
    ID = argv[7]
else:
    raise StandardError('Wrong number of args. Require 7')

# Make sure we have a valid constitutive model
if not( constit in ['vm', 'VM', 'H8', 'h8', 'anis', 'ANIS']):
    raise ValueError("Bad constit given '{}'.\nMust be 'vm', 'VM', 'H8', 'anis', 'ANIS'.".format(constit))

key = n.genfromtxt('ExptSummary.dat', delimiter=', ')
key = key[ key[:,0] == expt ].ravel()
a_true, R, t, X, sigma, torque, dmax, phimax, Lg = key[4:]
torque*=(2*pi*R*R*t) # torque to force
if constit in ['vm', 'VM']:
    # The *500 is b/c abaqus behaves odd when the cloads are O(1) for vm
    # The behavior is normal when the cloads are O(1000)
    # We'll divide the max LPF by 500 further down
    # I found that with H8 it is normal with cloads O(1)
    torque*=500
K = a_true/R # load ratio
force = K*torque # force
# Disp. control:  Monitor axial disp
# Can't monitor rotation (UR3 of a node is not rot'n about z-axis)
riks_DOF_num = 3
riks_DOF_val = 1.1*(dmax)/2

print('Cload 3 magnitude: {:.2f}'.format(force))
print('Cload 6 magnitude: {:.2f}'.format(torque))
print('Max displacement of Riks Node = {:.8f}'.format(riks_DOF_val))

# Load up the node and element lists
nodelist = n.load('./ConstructionFiles/abaqus_nodes.npy')
elemlist = n.load('./ConstructionFiles/abaqus_elements.npy')

fid =  open('{}.inp'.format(inpname),'w')

## Info for me
fid.write('** TTGM={}\n'.format(expt))
fid.write('** alpha = {:.3f}\n'.format(a_true))
fid.write('** num_el_fine_th = {}\n'.format(num_el_fine_th))
fid.write('** Imperfection = {} #Not percentage\n'.format(dt))
fid.write('** Eccentricity = {} #Not percentage\n'.format(ecc))
fid.write('** Rm = {:.4f}\n'.format(R))
fid.write('** tg = {:.4f}\n'.format(t))
fid.write('** ID = {}\n'.format(ID))
fid.write('** constit = "{}"\n'.format(constit))
z = time.localtime()
fid.write('** Generated on {}/{}/{} at {}:{}\n'.format(*z))

## HEADING and PREPRINT
fid.write('*Heading\n' +
          'Z is Tube Axis\n'  +
          'Min Wall Thickness is (x>0,y=0)\n'
          )
fid.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')

###################
## PART and sets
###################
fid.write('****************************************\n')
fid.write('***************** PART *****************\n')
fid.write('****************************************\n')
fid.write('*part, name=PART\n')
# Nodes
fid.write('*node, nset=NS_ALLNODES\n')
for i in nodelist:
    fid.write('{:.0f}, {:.12f}, {:.12f}, {:.12f}\n'.format(i[0],i[1],i[2],i[3]))
# Elements
fid.write('*element, type=C3D8R, elset=ES_ALLELEMENTS\n')
for i in elemlist:
    fid.write('{}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(
               i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8])
              )
# Node and Element Sets
with open('./ConstructionFiles/abaqus_sets.txt','r') as setfid:
    sets = setfid.read()
    fid.write(sets)
    setfid.close()

# One more dip-rot tracking node, that nearest to X=1.9685/2, Y=0, Z = 0.64
# I'm doing it here b/c it's easy to work with the entire nodelist
po = n.array([1.9675/2, 0, Lg/2])
loc = n.argmin( n.linalg.norm(nodelist[:,1:] - po, axis=1) )
nodenum = nodelist[loc, 0]
fid.write('*nset, nset=NS_DISPROT_NEW\n')
fid.write('{:.0f}\n'.format(nodenum))

# Orientation, transformation, section
fid.write('*orientation, name=ANISOTROPY, system=cylindrical, definition=coordinates\n' +
          '0, 0, 0, 0, 0, 1\n'
          )
fid.write('*solid section, elset=ES_ALLELEMENTS, material=MATERIAL, orientation=ANISOTROPY\n')
fid.write('*Hourglass Stiffness\n' +
          '40.0, , , \n'
          )
fid.write('*transform, nset=NS_ALLNODES, type=C\n' +
          '0, 0, 0, 0, 0, 1\n'
          )
fid.write('*end part\n')
# end part

###################
#### ASSEMBLY ####
###################
fid.write('****************************************\n')
fid.write('*************** Assembly ***************\n')
fid.write('****************************************\n')
fid.write('*assembly, name=ASSEMBLY\n')
fid.write('*instance, name=INSTANCE, part=PART\n')
fid.write('*end instance\n')
# Reference points
nodenum, zcoord = nodelist[:,0].max()+1, nodelist[:,3].max()
fid.write('*node, nset=NS_RPTOP\n' +
          '{:.0f}, 0, 0, {:.12f}\n'.format(nodenum, zcoord)
          )
# transform for RPs
fid.write('*transform, nset=NS_RPTOP, type=C\n' +
          '0, 0, 0, 0, 0, 1\n')
nodenum+=1
fid.write('*node, nset=NS_RPBOT\n' +
          '{:.0f}, 0, 0, 0\n'.format(nodenum)
          )
# Riks monitoring point, must be defined in the assembly
fid.write('** Riks displacement monitoring node must be defined in the assemlby\n' + 
          '*nset, nset=RIKSMON, instance=INSTANCE\n' + 
          'NS_DISPROT_NEW\n')
# Surfaces
fid.write('*surface, type=node, name=SURF_TOPSURFACE\n' +
          'INSTANCE.NS_TOPSURFACE\n'
          )
fid.write('*surface, type=node, name=SURF_BOTSURFACE\n' +
          'INSTANCE.NS_BOTTOMSURFACE\n'
          )
# Kinematic coupling
fid.write('*orientation, name=ORI_COUPLING, system=cylindrical, definition=coordinates\n' +
          '0, 0, 0, 0, 0, 1\n'
          )
fid.write('*coupling, constraint name=CONSTRAINT_TOPSURFACE, ref nod=NS_RPTOP, surface=SURF_TOPSURFACE, orientation=ORI_COUPLING\n' +
          '*kinematic\n'
          )
fid.write('*coupling, constraint name=CONSTRAINT_BOTSURFACE, ref nod=NS_RPBOT, surface=SURF_BOTSURFACE, orientation=ORI_COUPLING\n' +
          '*kinematic\n' +
          '2, 3\n' +
          '5, 6\n'
          )
fid.write('*end assembly\n')

# end assembly

###################
#### MATERIAL #####
###################
fid.write('****************************************\n')
fid.write('***************  MATERIAL **************\n')
fid.write('****************************************\n')
fid.write('*material, name=MATERIAL\n')
if constit in ['vm', 'VM']:
    with open('./ConstructionFiles/abaqus_material_VM_TT20.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()
elif constit in ['h8','H8']:
    with open('./ConstructionFiles/abaqus_material_H8_TT20.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()    
elif constit in ['ANIS', 'anis']:
    with open('./ConstructionFiles/abaqus_material_ANIS.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()

###################
### INITIAL BCs ###
###################
fid.write('****************************************\n')
fid.write('**************  INITIAL BCs ************\n')
fid.write('****************************************\n')
fid.write('*boundary\n' +
          'ASSEMBLY.NS_RPBOT, 1, 6 \n' + 
          'ASSEMBLY.NS_RPTOP, 1, 2 \n' +      # Top ref point can only translate up and rotate about axis
          'ASSEMBLY.NS_RPTOP, 4, 5 \n'
          )

###################
###### STEP #######
###################
fid.write('****************************************\n')
fid.write('****************** STEP ****************\n')
fid.write('****************************************\n')
fid.write('*step, name=STEP, nlgeom=yes, inc=500\n')
# [1]Inital arc len, [2]total step, [3]minimum increm, [4]max increm (no max if blank), [5]Max LPF, [6]Node whose disp is monitored, [7]DOF, [8]Max Disp
if not n.isnan(a_true):
    # Riks if tension and torsion
    LPFmax = 1.1
    if constit in ['vm','VM']:
        LPFmax*=(1/500)
        fid.write('*static, riks\n' +
                '0.001, 1.0, 1e-10, .005, {:.4f}, ASSEMBLY.RIKSMON, {:.0f}, {:.6f}\n'.format(LPFmax,riks_DOF_num, riks_DOF_val)
                  )
    else:
        fid.write('*static, riks\n' +
                '0.075, 1.0, 1e-10, .075, {:.4f}, ASSEMBLY.RIKSMON, {:.0f}, {:.6f}\n'.format(LPFmax,riks_DOF_num, riks_DOF_val)
                  )
    fid.write('**[1]Inital arc len, [2]total step, [3]minimum increm, [4]max increm (no max if blank), [5]Max LPF, [6]Node whose disp is monitored, [7]DOF, [8]Max Disp\n')
    fid.write('*cload\n' +
              'ASSEMBLY.NS_RPTOP, 3, {:.5f}\n'.format(force) + 
              'ASSEMBLY.NS_RPTOP, 6, {:.5f}\n'.format(torque)
              )
    fid.write('** Increase number of permitted cutbacks in increment\n'
              '*controls, parameters=time incrementation\n' + 
              ' , , , , , , , 10, , , , , \n'
              )
else:
    fid.write('*static\n' +
              '0.005, 1., 1e-05, .005\n'
              )
              

# field output
fid.write('*output, field, frequency=1\n')
fid.write('** COORn must be called under history output, but COORD can be called in field\n')
fid.write('*node output, nset=INSTANCE.NS_DISPROT_LO\n' +   # disprot nodesets
          'U, UR, COORD\n'
          )
fid.write('*node output, nset=INSTANCE.NS_DISPROT_HI\n' +
          'U, UR, COORD\n'
          )
fid.write('*node output, nset=INSTANCE.NS_DISPROT_NEW\n' +
          'U, UR, COORD\n'
          )
fid.write('*node output, nset=ASSEMBLY.NS_RPTOP\n' +    # refpt node
          'U, UR, CF\n'
          )
fid.write('*node output, nset=ASSEMBLY.NS_RPBOT\n' +
          'RF, RM\n'
          )          
fid.write('*node output, nset=INSTANCE.NS_RADIALCONTRACTION\n' +    # radial contraction set
          'U, UR\n'
          )
for i in ['ES_Z', 'ES_THICKNESS', 'ES_THICKNESS_BACK', 'ES_ANALZONE']:
    fid.write('*element output, elset=INSTANCE.{}, directions=YES\n'.format(i) +    # sts, stn in element sets
              'S, PE, LE, COORD'
              )
    if constit in ['H8','h8','anis','ANIS']:
        fid.write(', SDV1, SDV2\n')
    else:
        fid.write(', PEEQ, MISESONLY\n')

fid.write('*end step\n')
# end step

'''
Abaqus Users Guide: Sxn 2.1.5
Output database output of field vector-valued quantities at transformed nodes is in the global system. The local transformations are also written to the output database. You can apply these transformations to the results in the Visualization module of Abaqus/CAE to view the vector components in the transformed systems.
'''
fid.close()


'''
# History output (coor)
fid.write('** COORn must be called under history output\n')
fid.write('*output, history, frequency=1\n')
fid.write('*node output, nset=INSTANCE.NS_DISPROT_LO\n' +
          'COOR1, COOR2, COOR3\n'
          )
fid.write('*node output, nset=INSTANCE.NS_DISPROT_HI\n' +
          'COOR1, COOR2, COOR3\n'
          )
fid.write('*node output, nset=INSTANCE.NS_RADIALCONTRACTION\n' +    # radial contraction set
          'COOR1, COOR2, COOR3\n'
          )        
'''
