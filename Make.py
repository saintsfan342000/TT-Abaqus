import numpy as n
import os

try:
    argv, expt, num_el_fine_th, dt, eccen = argv
    d = n.genfromtxt('ExptParameters.dat', delimiter=', ', dtype=str)
    # [0]Expt no, [1]IDg, [2]Mean thickness, [3]Min thickness
    ID, tg, tmin = d[ d[:,0]==expt, 1:].ravel()   
    if eccen = 'auto':
        eccen = 1 - tmin/tg
except:
    ID = '1.75' # Inner Radius
    tg = '0.038'
    num_el_fine_th = '12' # Num elements thru the test section thicknes
    dt = '.01'

Lg = '0.2'
Ltop = '1.3'  # Length of thick section above radius/chamf
ODtop = '1.9685'    # Radius of thick section    
R = '0.125'    # Radius of chamf
    
### Nodes
os.system('python Nodes.py {:s} {:s} {:s} {:s} {:.4f} {.3f} {:d} {.3f}'.format(
                Lg, Ltop, ODtop, ID, tg, R, num_el_fine_th, dt)

### Elements
os.system('python Elements.py')

### Sets
os.system('python Sets.py')

### Eccentricity

