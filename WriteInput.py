import numpy as n
from numpy import (linspace, asarray, in1d, empty,
                    hstack , vstack, pi, compress)
from sys import argv

'''
Writes the input file
'''

try:
    inpname = argv[1]
    alpha = float(argv[2])
    force = float(argv[3])
    torque = float(argv[4])
    riks_DOF_num = int(argv[5])
    riks_DOF_val = float(argv[6])
    constit = argv[7]
except IndexError:
    inpname = 'Input2'
    alpha = 0.5
    force = 1000
    torque = 2000
    riks_DOF_num = 6
    riks_DOF_val = .05
    constit = 'vm'

nodelist = n.load('./ConstructionFiles/abaqus_nodes.npy')
elemlist = n.load('./ConstructionFiles/abaqus_elements.npy')

fid =  open('{}.inp'.format(inpname),'w')

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
nodenum, zcoord = nodelist[:,0].max()+1, nodelist[:,2].max()
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
          'NS_DISPROT_LO\n')
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

if constit in ['vm', 'VM']:
    with open('./ConstructionFiles/abaqus_material_VM.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()
elif constit == 'H8':
    with open('./ConstructionFiles/abaqus_material_H8.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()    
elif constit in ['ANIS', 'anis']:
    with open('./ConstructionFiles/abaqus_material_ANIS.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()
else:
    raise ValueError("Bad constit given '{}'.\nMust be 'vm', 'VM', 'H8', 'anis', 'ANIS'.".format(constit))

###################
### INITIAL BCs ###
###################
fid.write('*boundary\n' +
          'ASSEMBLY.NS_RPBOT, 1, 6 \n' + 
          'ASSEMBLY.NS_RPTOP, 1, 2 \n' +      # Top ref point can only translate up and rotate about axis
          'ASSEMBLY.NS_RPTOP, 4, 5 \n'
          )

###################
###### STEP #######
###################

fid.write('*step, name=STEP, nlgeom=yes, inc=500\n')
# [1]Inital arc len, [2]total step, [3]minimum increm, [4]max increm (no max if blank), [5]Max LPF, [6]Node whose disp is monitored, [7]DOF, [8]Max Disp
if alpha != 'PS':
    # Riks if tension and torsion
    fid.write('*static, riks\n' +
            '0.01, 1.0, 1e-05, .005, 2, ASSEMBLY.RIKSMON, {:.0f}, {:.3f}\n'.format(riks_DOF_num, riks_DOF_val)
              )
    fid.write('**[1]Inital arc len, [2]total step, [3]minimum increm, [4]max increm (no max if blank), [5]Max LPF, [6]Node whose disp is monitored, [7]DOF, [8]Max Disp\n')
    fid.write('*cload\n' +
              'ASSEMBLY.NS_RPTOP, 3, {:.2f}\n'.format(force) + 
              'ASSEMBLY.NS_RPTOP, 6, {:.2f}\n'.format(torque)
              )
else:
    fid.write('*static\n' +
              '0.005, 1., 1e-05, .005\n'
              )
              

# field output
fid.write('*output, field, frequency=1\n')
fid.write('*node output, nset=INSTANCE.NS_DISPROT_LO\n' +   # disprot nodesets
          'U, UR\n'
          )
fid.write('*node output, nset=INSTANCE.NS_DISPROT_HI\n' +
          'U, UR\n'
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
for i in ['ES_Z', 'ES_TH', 'ES_TH_BACK']:
    fid.write('*element output, elset=INSTANCE.{}, directions=YES\n'.format(i) +    # sts, stn in element sets
              'S, PE, LE'
              )
    if constit in ['H8','anis','ANIS']:
        fid.write(', SDV1, SDV2\n')
    else:
        fid.write('\n')
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
fid.write('*end step\n')
# end step

'''
Abaqus Users Guide: Sxn 2.1.5
Output database output of field vector-valued quantities at transformed nodes is in the global system. The local transformations are also written to the output database. You can apply these transformations to the results in the Visualization module of Abaqus/CAE to view the vector components in the transformed systems.
'''
fid.close()
