## usage: /Users/lixinmeng/anaconda2/bin/python *.py
## notice : there are two missing coupling points
## !!!!! array, coloum  3 31
## !!!!! array, coloum  22 14
## ####################################

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib as mplt

InFile = 'inner-inter-state.txt' #'states-angle-2d-noborder.txt'


mix = []
dif = []


with open(InFile,'r') as f:
    rawdata = f.readlines()
    for line in rawdata[:]:

        list = line.split()
        ## list[0], list[1], list[2] split it
        
        #indexi = int(list[0]) -1 
        #indexj = int(list[1]) -1
        #index1 = (indexi % 56)
        #index2 = (indexj % 56)
        
        #grida = int( (index1+index2)/2.0 )
        #gridb = int( max( indexi  ,indexj ) / 56.0 )
        t = float(list[2])
        mix.append(t)
        #result[frame][gridb][grida] = t

# 78400 /2 
for j in xrange(39200):
    demo = mix[j*2] - mix[j*2+1]
    dif.append(demo)


### plot distribution


fig5, ax5 = plt.subplots()


#ax5.set_xlim(-70.0,70.0)
#ax5.set_xticks(np.linspace(-60,60,7,endpoint=True))
#ax5.set_ylim(0,0.04)
#ax5.set_yticks(np.linspace(0,0.04,4,endpoint=False))
#fig5.set_size_inches(6, 3)

#### http://xkcd.com/color/rgb/

sns.distplot( dif,bins=30,norm_hist=True, kde = False, color = "red",label = 'neighboring-inter-intra-difference'   )
#sns.distplot( inner,bins=500 ,  norm_hist= True, kde = False, color = sns.xkcd_rgb["pale red"] ,label = 'intra dimer'   )


#####
# the result here is 
# mu=-21.0880697245, sigma=11.6533229858
# mu=23.5993130026, sigma=13.3824014869 
####

#plt.xlim(0.65,0.85)
#plt.xticks ( [0.65, 0.739, 0.85])
#plt.xticks(np.linspace(0.65,0.85,5,endpoint=True))

#plt.ylim(0,0.04)
#plt.yticks(np.linspace(0,0.04,0,endpoint=True))

plt.xlabel('Angle ($^\circ$)')
plt.ylabel('Probability Density')

plt.legend(loc='best')


#plt.legend(loc='best')
fig5.savefig('minimal-model-inner-inter-dif.pdf',bbox_inches='tight')

#
#mu=19.9064230234, sigma=14.5756637612
#mu=-22.9004333214, sigma=16.4392568036
