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

readstsstn=True

try:
    job = argv[1].split('.')[0] # Job name, cutting off the .odb extension in case it was given
except IndexError:
    job = 'Input3'
    from glob import glob
    job = glob('*.odb')[0].split('.')[0]

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

'''
# Transformation of field values
# Apparently must be done for each field value for each frame
csysname = h_odb.rootAssembly.datumCsyses.keys()[0] #'ASSEMBLY_INSTANCE_ANISOTROPY'
h_csys = h_odb.rootAssembly.datumCsyses[csysname]
fv = frame.fieldOutputs['fv_key']   # No .values!
transformed_fv = fv.getTransformedField(
    datumCsys=h_csys)
'''

h_nset_rp_top = h_odb.rootAssembly.nodeSets[nset_rp_top]
h_nset_dr_lo = h_inst.nodeSets[nset_dr_lo]
h_nset_dr_hi = h_inst.nodeSets[nset_dr_hi]
h_nset_dr_new = h_inst.nodeSets[nset_dr_new]
h_elset_th = h_inst.elementSets[elset_th]

F = np.empty( (num_incs) )
T = np.empty( (num_incs) )
d_lo = np.empty( (num_incs, 3) )
d_hi = np.empty( (num_incs, 3) )
d_new = np.empty( (num_incs, 3) )

# Grab undef coords of dr_lo and hi
c_lo = h_nset_dr_lo.nodes[0].coordinates[:]
Lg_lo = c_lo[2]
c_hi = h_nset_dr_hi.nodes[0].coordinates[:]
Lg_hi = c_hi[2]
c_new = h_nset_dr_new.nodes[0].coordinates[:]
Lg_new = c_new[2]

if readstsstn:
    S = np.empty((num_incs,6))
    LE = np.empty_like(S)
    PE = np.empty_like(S)
    # For the local coordinate system of the node in center of ES_TH
    LCS = np.empty((num_incs,3,3))

for i in range(num_incs):
    frame = h_All_Frames[i]
    F[i] = frame.fieldOutputs['CF'].getSubset(region=h_nset_rp_top).values[0].data[2]
    T[i] = frame.fieldOutputs['CM'].getSubset(region=h_nset_rp_top).values[0].data[2]
    d_lo[i] = frame.fieldOutputs['U'].getSubset(region=h_nset_dr_lo).values[0].data[:]
    d_hi[i] = frame.fieldOutputs['U'].getSubset(region=h_nset_dr_hi).values[0].data[:]
    d_new[i] = frame.fieldOutputs['U'].getSubset(region=h_nset_dr_new).values[0].data[:]
    if readstsstn:
        # Get S for each element in the elset and take the mean
        fv = frame.fieldOutputs['S'].getSubset(region=h_elset_th).values
        S[i] = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)
        for j in range(3):
            LCS[i,:,j] = fv[4].localCoordSystem[j]
        # Log stn
        fv = frame.fieldOutputs['LE'].getSubset(region=h_elset_th).values
        LE[i] = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)
        # Plastic stn
        fv = frame.fieldOutputs['PE'].getSubset(region=h_elset_th).values
        PE[i] = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)

h_odb.close()

# Calc phi and delta
for i in ['_lo','_hi','_new']:
    # Add on undef coords
    exec("d%s+=c%s"%(i,i))
    # phi
    exec("phi%s = 2*(180/pi)*np.arctan2(d%s[:,1],d%s[:,0])"%(i,i,i))
    exec("phi%s-=phi%s[0]"%(i,i))
    # Delta
    exec("d%s = (d%s[:,2]-d%s[0,2])/d%s[0,2]"%(i,i,i,i))

def headerline(fname, hl):
    fid = open(fname, 'r')
    data = fid.read()
    fid.close()
    del fid
    fid = open(fname, 'w')
    if hl[0] != '#':
        fid.write('#')
    fid.write(hl)
    if hl[-1] != '\n':
        fid.write('\n')
    fid.write(data)
    fid.close()

# Save
np.savetxt('disprot_new.dat',X = np.vstack((F, T, d_lo, phi_lo, d_hi, phi_hi, d_new, phi_new)).T, fmt='%.12f', delimiter=', ')
headerline('disprot_new.dat','#[0]Nom AxSts, [1]Nom Shear Sts, [2]d/Lg lo, [3]Rot lo, [4]d/Lg hi, [5]Rot hi, [6]d/Lg new, [7]Rot new\n')

if readstsstn:
    np.savetxt('S.dat', X=S, fmt='%.6f', delimiter=',')
    headerline('S.dat','[0]Srr, [1]Sqq, [2]Szz, [4]Srq, [5]Srz?, [6]Sqz')
    np.savetxt('LE.dat', X=LE, fmt='%.6f', delimiter=',')
    headerline('LE.dat','[0]LErr, [1]LEqq, [2]LEzz, [4]LErq, [5]LErz?, [6]LEqz')
    np.savetxt('PE.dat', X=PE, fmt='%.6f', delimiter=',')
    headerline('PE.dat','[0]PErr, [1]PEqq, [2]PEzz, [4]PErq, [5]PErz?, [6]PEqz')
    np.save('LCS.npy', LCS)
