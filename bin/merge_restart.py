# AUTHOR        : Caglar
# DATE          : Jan30,2020
import abfile
import numpy as np
ABphy = '/cluster/work/users/cagyum/TP5a0.06/expt_01.1/data/restart.2015_213_00_0000_copy.a'
ABeco = '/cluster/work/users/cagyum/TP5a0.06/expt_01.1/data/restart.2018_186_00_0000.a'

first_eco_variable = 'CO2_dic' # inspect the eco restart file and assign the first ecomodel variable
I1 = 760; I2 = 800
rABphy = abfile.ABFileRestart(ABphy,"r",idm=I2,jdm=I1)
rABeco = abfile.ABFileRestart(ABeco,"r",idm=I2,jdm=I1)

newfile = 'merged'
new_abfile = abfile.ABFileRestart(newfile,"w",idm=I2,jdm=I1)
new_abfile.write_header(rABphy._iexpt,rABphy._iversn,rABphy._yrflag,rABphy._sigver,rABphy._nstep,rABphy._dtime,rABphy._thbase)

for phykeys in sorted( rABphy.fields.keys() ) :
    fieldname = rABphy.fields[phykeys]["field"]
    k         = rABphy.fields[phykeys]["k"]
    t         = rABphy.fields[phykeys]["tlevel"]
    field     = rABphy.read_field(fieldname,k,t)
    print("Copying %10s to level %3d at time=%d"%(fieldname,k,t))
    new_abfile.write_field(field,None,fieldname,k,t)

for ecokeys in sorted( rABeco.fields.keys() ) :
    fieldname = rABeco.fields[ecokeys]["field"]
    k         = rABeco.fields[ecokeys]["k"]
    t         = rABeco.fields[ecokeys]["tlevel"]
    if fieldname == first_eco_variable:
       if (ecokeys == phykeys + 1):
           execute = 'True'
       else:
           print('expected ecosystem model key does not match the source ecosystem restart file key')
           print('this is possibly due to different number of variables in source physics and ecosystem physics')
           
           answer = input('Do you still want to continue? "y" or "n"  :')
           if answer == 'y' or answer == 'Y' or answer == 'yes' or answer == 'YES':
              execute = 'True'
              break
           else:
              execute = 'False'
              break
       break

if (execute):
   for keys in range(ecokeys,np.max(rABeco.fields.keys())+1):
       fieldname = rABeco.fields[keys]["field"]
       k         = rABeco.fields[keys]["k"]
       t         = rABeco.fields[keys]["tlevel"]
       field     = rABeco.read_field(fieldname,k,t)
       print("Copying %10s to level %3d at time=%d"%(fieldname,k,t))
       new_abfile.write_field(field,None,fieldname,k,t)

new_abfile.close()
