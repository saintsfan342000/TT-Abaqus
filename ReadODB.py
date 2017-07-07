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
from glob import glob
pi = np.pi

readstsstn = True
Anal_Zone = False

try:
    job = argv[1].split('.')[0] # Job name, cutting off the .odb extension in case it was given
except IndexError:
    odbs = glob('*.odb')
    if len(odbs)>1:
        raise ValueError('More than one odb in this directory, cant glob the filename.')
    job = odbs[0].split('.')[0]

inputs = glob('*.inp')
if len(inputs) != 1:
    print '%d .inp files in pwd; can only glob if theres 1.  Assuming default values for Rm and tg'%len(inputs)
    Rm = .8338
    tg = 0.0389
else:
    fid = open(inputs[0], 'r')
    count = 0
    for k,line in enumerate(fid):
        if line.rfind('Rm ') != -1:
            Rm = float( line.split(' ')[-1] )
            count += 1
        if line.rfind('tg =') != -1:
            tg = float( line.split(' ')[-1] )
            count += 1
        if count == 2:
            break

if not os.path.isfile(job + '.odb'):
    raise ValueError('The specified job name "%s" does not exist.'%(job))

nset_rp_top = 'NS_RPTOP'
nset_rp_bot = 'NS_RPBOT'
nset_dr_lo = 'NS_DISPROT_LO'
nset_dr_hi = 'NS_DISPROT_HI'
nset_dr_new = 'NS_DISPROT_NEW'
elset_th = 'ES_THICKNESS'
elset_analzone = 'ES_ANALZONE'

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
fv = h_ffos['fv_key']   # No .values!
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
    csysname = h_odb.rootAssembly.datumCsyses.keys()[0]
    h_csys = h_odb.rootAssembly.datumCsyses[csysname]
    S = np.empty((num_incs,6))
    LE = np.empty_like(S)
    PE = np.empty_like(S)
if Anal_Zone:
    csysname = h_odb.rootAssembly.datumCsyses.keys()[0]  
    h_csys = h_odb.rootAssembly.datumCsyses[csysname]
    h_elset_analzone = h_inst.elementSets[elset_analzone]
    S_Anal_Zone = np.empty((num_incs,6))

for i in range(num_incs):
    h_ffos = h_All_Frames[i].fieldOutputs  # A handle for all this frame's field outputs
    F[i] = h_ffos['CF'].getSubset(region=h_nset_rp_top).values[0].data[2]
    T[i] = h_ffos['CM'].getSubset(region=h_nset_rp_top).values[0].data[2]
    d_lo[i] = h_ffos['U'].getSubset(region=h_nset_dr_lo).values[0].data[:]
    d_hi[i] = h_ffos['U'].getSubset(region=h_nset_dr_hi).values[0].data[:]
    d_new[i] = h_ffos['U'].getSubset(region=h_nset_dr_new).values[0].data[:]
    if readstsstn:
        # Get S for each element in the elset and take the mean
        fv = h_ffos['S'].getTransformedField(datumCsys=h_csys).getSubset(region=h_elset_th).values
        S[i] = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)
        # Log stn
        fv = h_ffos['LE'].getTransformedField(datumCsys=h_csys).getSubset(region=h_elset_th).values
        LE[i] = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)
        # Plastic stn
        fv = h_ffos['PE'].getTransformedField(datumCsys=h_csys).getSubset(region=h_elset_th).values
        PE[i] = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)
    if Anal_Zone:
        # Get S for each element in the elset and take the mean
        fv = h_ffos['S'].getgetSubset(region=h_elset_analzone).getSubset(region=h_elset_analzone).values
        S_Anal_Zone = np.array([ fv[k].data for k in range(len(fv)) ]).mean(axis=0)

        
h_odb.close()

# F and T to nom sts
F *= (1/(2*pi*Rm*tg))
T *= (1/(2*pi*Rm*Rm*tg))

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
np.savetxt('disprot.dat',X = np.vstack((F, T, d_lo, phi_lo, d_hi, phi_hi, d_new, phi_new)).T, fmt='%.12f', delimiter=', ')
headerline('disprot.dat','#[0]Nom AxSts, [1]Nom Shear Sts, [2]d/Lg lo, [3]Rot lo, [4]d/Lg hi, [5]Rot hi, [6]d/Lg new, [7]Rot new\n')

if readstsstn:
    np.savetxt('S.dat', X=S, fmt='%.6f', delimiter=',')
    headerline('S.dat','[0]Srr, [1]Sqq, [2]Szz, [4]Srq, [5]Srz?, [6]Sqz')
    np.savetxt('LE.dat', X=LE, fmt='%.6f', delimiter=',')
    headerline('LE.dat','[0]LErr, [1]LEqq, [2]LEzz, [4]LErq, [5]LErz?, [6]LEqz')
    np.savetxt('PE.dat', X=PE, fmt='%.6f', delimiter=',')
    headerline('PE.dat','[0]PErr, [1]PEqq, [2]PEzz, [4]PErq, [5]PErz?, [6]PEqz')
