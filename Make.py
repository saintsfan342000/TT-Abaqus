import numpy as n
pi = n.pi
import os
from sys import argv, exit

def println():
    print('-----------------------------------------------------------')

def printuse():
    println()
    print('Requires 6 command-line arguments.')
    print('[1] Experiment no. whose geometry and loading you\'re modeling.')
    print('[2] Number of elements thru-thickness in the test section')
    print('[3] Imperfection magnitude (as fraction, not percentage)')
    print('[4] Eccentricity.  Give a number, or give "auto" and the geometry will match that of the test specimen.')
    print('[5] Input file name')
    print('[6] Constitutive model: "vm", "h8", or "anis"')
    println()


if (len(argv) != 7) and (len(argv) != 2):
    printuse()
    exit()
elif (len(argv) == 2):
    if argv[1] in ['TestVM', 'TestH8']:
        (expt, num_el_fine_th, dt,
                eccen, inpname, constit) = ('2', '6', '0',
                                            'auto', argv[1], argv[1][-2:])
    else:
        printuse()
        exit()
else:
    argv, expt, num_el_fine_th, dt, eccen, inpname, constit = argv


d = n.genfromtxt('ExptSummary.dat', delimiter=',', dtype=str)
# [0]Expt no, [1]IDg, [2]Mean thickness, [3]Min thickness
alpha, Rm, tg, *X = d[ d[:,0]==expt, 4:].ravel()   
X = X[0]  # Because now I append extra info to ExptSumm
ID = '{:.4f}'.format((float(Rm)-float(tg)/2)*2)
if eccen == 'auto':
    eccen = float(X)/100
    print('Eccen = {:.2f}'.format(eccen*100))

inpname = inpname.split('.')[0] # In case I give an extension.

Lg = '0.2'
Ltop = '.5'  # Length of thick section above radius/chamf
ODtop = '1.9675'    # Radius of thick section    
R = '0.125'    # Radius of chamf

### Nodes
println()
print('Generating nodes.')
os.system('python Nodes.py {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(
            Lg, Ltop, ODtop, ID, tg, R, num_el_fine_th, dt))
print('Done.')

### Elements
println()
print('Generating elements.')
os.system('python Elements.py')
print('Done.')

### Sets
println()
print('Generating node and element sets.')
os.system('python Sets.py')
print('Done.')

### Eccentricity
println()
print('Generating Eccentricity.')
os.system('python Eccen.py {}'.format(eccen))
print('Done.')

### Write Input
println()
print('Writing input file, {}.inp'.format(inpname))
os.system('python WriteInput.py {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(
          expt, inpname, constit, num_el_fine_th, 
          dt, str(eccen), ID ))
print('Done!')
