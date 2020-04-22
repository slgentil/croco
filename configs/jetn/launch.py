import crocosi.launch as lc
import numpy as np

#lc.run('test', launch=True)

# 1 run of 100 days
#lc.run('fsturb_4km_0a100j'.format(m), launch=True)

# multiple runs, here multiple modes
R = {0: {'fsturb_avmodes': [1., 0., 0., 0., 0.]},
     1: {'fsturb_avmodes': [0., 1., 0., 0., 0.]}
    }
#for m, mp in R.items(): 
#    lc.run('fsturb_m{}_4km_0a1000j'.format(m), jobname='{}'.format(m),
#            nbchains=10, launch=True, **mp)

# change amplitude
#fsturb: dt     F     Lmin     Lmax   t_start  t_decay  Nmode xmid xwid   ymid     ywid
#        10.  1.e-3  50.e+3   100.e+3      0.      800.     5    0.    0. 1440.e+3  500.e+3
# should be a dict
def get_fsturb(f):
    return [10., f, 50.e+3, 100.e+3,0., 800., 
             5, 0., 0., 1440.e+3, 500.e+3]
F = np.arange(1.e-4, 2.e-3, 5.e-4)
for m, mp in R.items():
    for fi, f in enumerate(F):
        lc.run('fsturb_m{}_a{}_4km_0a300j'.format(m,fi), 
               jobname='{}f{}'.format(m,fi),
               nbchains=3, launch=True, fsturb=get_fsturb(f),**mp)



