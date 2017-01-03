import numpy as n
pi = n.pi
import os
from sys import argv

try:
    argv, expt, num_el_fine_th, dt, eccen, inpname, constit = argv
    d = n.genfromtxt('ExptParams.dat', delimiter=', ', dtype=str)
    # [0]Expt no, [1]IDg, [2]Mean thickness, [3]Min thickness
    alpha, ID, tg, tmin = d[ d[:,0]==expt, 1:].ravel()   
    if eccen == 'auto':
        eccen = 1 - float(tmin)/float(tg)
        print('Eccen = {:.2f}'.format(eccen*100))
except ValueError:
    raise ValueError('Problem!')
    expt = '20'
    ID = '1.75' # Inner Radius
    tg = '0.038'
    num_el_fine_th = '12' # Num elements thru the test section thicknes
    dt = '.01'
    inpname = 'Input_20'
    constit = 'vm'

Lg = '0.2'
Ltop = '1.3'  # Length of thick section above radius/chamf
ODtop = '1.9685'    # Radius of thick section    
R = '0.125'    # Radius of chamf
    
### Nodes
os.system('python Nodes.py {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(
                Lg, Ltop, ODtop, ID, tg, R, num_el_fine_th, dt))

### Elements
os.system('python Elements.py')

### Sets
os.system('python Sets.py')

### Eccentricity
os.system('python Eccen.py {}'.format(eccen))

### Write Input
os.system('python WriteInput.py {:s} {:s} {:s}'.format(expt, inpname, constit))
