''' Read ODB
job_name:   ODB name without odb extension
inst_name:  name of instance
step_name:  name of load step
ASSEMBLY.NS_RPTOP : The top reference point to which force and torque are applied
ASSEMBLY.NS_RPBOT : The bottom ref. pt. which is fixed
INSTANCE.NS_DISPROT_LO : The node just on the cusp of the the test section
INSTANCE.NS_DISPROT_HI : The node one above LO (in case I want to interpolate btwn them)
INSTANCE.ES_Z : The line of elements at theta=0
'''

import odbAccess as O
import numpy as np
import os
from sys import argv
pi = np.pi

try:
    job = argv[1].split('.')[0] # Job name, cutting off the .odb extension in case it was given
except IndexError:
    job = 'Input3'

if not os.path.isfile(job + '.odb'):
    raise ValueError('The specified job name "%s" does not exist.'%(job))

nset_rp_top = 'NS_RPTOP'
nset_rp_bot = 'NS_RPBOT'
nset_dr_lo = 'NS_DISPROT_LO'
nset_dr_hi = 'NS_DISPROT_HI'
nset_dr_new = 'NS_DISPROT_NEW'
elset_th = 'ES_THICKNESS'

h_odb = O.openOdb(job + '.odb',readOnly=True)
h_inst = h_odb.rootAssembly.instances[ h_odb.rootAssembly.instances.keys()[0] ]
h_step = h_odb.steps[ h_odb.steps.keys()[0] ]
h_All_Frames = h_step.frames
num_incs = len(h_All_Frames)

h_nset_rp_top = h_odb.rootAssembly.nodeSets[nset_rp_top]
h_nset_dr_lo = h_inst.nodeSets[nset_dr_lo]
h_nset_dr_hi = h_inst.nodeSets[nset_dr_hi]
h_nset_dr_new = h_inst.nodeSets[nset_dr_new]
h_elset_th = h_inst.elementSets[elset_th]

F = np.empty( (num_incs) )
T = np.empty( (num_incs) )
d_lo = np.empty( (num_incs) )
d_hi = np.empty( (num_incs) )

# Grab undef coords of dr_lo and hi
Lg_lo = h_nset_dr_lo.nodes[0].coordinates[2]
Lg_hi = h_nset_dr_hi.nodes[0].coordinates[2]
Lg_new = h_nset_dr_new.nodes[0].coordinates[2]

dr_lo = np.empty((num_incs,3))
dr_hi = np.empty((num_incs,3))
dr_new = np.empty((num_incs,3))

for i in range(num_incs):
    F[i] = h_All_Frames[i].fieldOutputs['CF'].getSubset(region=h_nset_rp_top).values[0].data[2]
    T[i] = h_All_Frames[i].fieldOutputs['CM'].getSubset(region=h_nset_rp_top).values[0].data[2]
    d_lo[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_lo).values[0].data[2]
    d_hi[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_hi).values[0].data[2]
    dr_lo[i] = h_All_Frames[i].fieldOutputs['COORD'].getSubset(region=h_nset_dr_lo).values[0].data[:]
    dr_hi[i] = h_All_Frames[i].fieldOutputs['COORD'].getSubset(region=h_nset_dr_hi).values[0].data[:]
    dr_new[i] = h_All_Frames[i].fieldOutputs['COORD'].getSubset(region=h_nset_dr_new).values[0].data[:]

h_odb.close()

phi_lo = 2*(180/pi)*np.arctan2(dr_lo[:,1],dr_lo[:,0])
phi_lo-=phi_lo[0]

phi_hi = 2*(180/pi)*np.arctan2(dr_hi[:,1],dr_hi[:,0])
phi_hi-=phi_hi[0]

phi_new = 2*(180/pi)*np.arctan2(dr_new[:,1],dr_new[:,0])
phi_new-=phi_new[0]

# Convert disp to d/Lg
d_lo, d_hi = d_lo/Lg_lo, d_hi/Lg_hi
d_new = (dr_new[:,2] - dr_new[0,2])/Lg_new

def headerline(fname, hl):
    fid = open(fname, 'r')
    data = fid.read()
    fid.close()
    del fid
    fid = open(fname, 'w')
    fid.write(hl)
    if hl[-1] != '\n':
        fid.write('\n')
    fid.write(data)
    fid.close()

# Save
np.savetxt('disprot.dat',X = np.vstack((F, T, d_lo, phi_lo, d_hi, phi_hi, d_new, phi_new)).T, fmt='%.12f', delimiter=', ')
headerline('disprot.dat','#[0]Nom AxSts, [1]Nom Shear Sts, [2]d/Lg lo, [3]Rot lo, [4]d/Lg hi, [5]Rot hi, [6]d/Lg new, [7]Rot new\n')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
