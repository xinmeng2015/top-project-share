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

InFile = 'states-angle-2d-no-border.txt' #'states-angle-2d-noborder.txt'



#result = np.zeros((28, 55))
inner = []
inter = []


NFRAMES = 50
LINES = 1566
## one frame 1566 lines
## total 78300/1566 = 50 
## frames = 50  i in frames:[0,49]  ==> 500 +(i+1)*10  ps
## 

with open(InFile,'r') as f:
    rawdata = f.readlines()
    for t in xrange(NFRAMES):
        time = 500 +(t+1)*10
        print time

        result = np.zeros((28, 55))
        for i in xrange(28):
            for j in xrange(55):
                #result[i][j] = 600.0
                result[i][j] = np.nan

        for line in rawdata[t*LINES:(t+1)*LINES]:
            #print (line)
            list = line.split()
            ## reorgnize indices

            ori_i = int(list[0])
            ori_j = int(list[1])

            indexi = ori_i - 60 - (int(ori_i/60)*4+2)
            indexj = ori_j - 60 - (int(ori_j/60)*4+2)

            indexi = indexi -1
            indexj = indexi -1

            index1 = (indexi % 56)
            index2 = (indexj % 56)
            grida = int( (index1+index2)/2.0 )
            gridb = int( max( indexi  ,indexj ) / 56.0 )


            t = float(list[2])
            result[gridb][grida] = t

        data = pd.DataFrame( result )
        fig, ax = plt.subplots()

        #sns.set_context("paper", font_scale=2.0)
        
        sns.heatmap(data,xticklabels=False, yticklabels=False,vmin=-60.0, vmax =60.0 )
        fig.set_size_inches(9, 3)
        #plt.xlabel(str(time)+'  ps')
        fig.suptitle(str(time)+'  ps', fontsize=20)

        file_name = str(time)+'-frame-heatmap.png'
        fig.savefig(file_name,bbox_inches='tight', dpi=360)



"""
        list = line.split()
        ## list[0], list[1], list[2] split it
        indexi = int(list[0]) -1 
        indexj = int(list[1]) -1
        t = float(list[2])
        diff = np.absolute(indexi-indexj)
        if (diff == 1):
            inner.append(t)
        else:
            inter.append(t)     
 


### plot distribution


fig5, ax5 = plt.subplots()


#ax5.set_xlim(-70.0,70.0)
#ax5.set_xticks(np.linspace(-60,60,7,endpoint=True))
#ax5.set_ylim(0,0.04)
#ax5.set_yticks(np.linspace(0,0.04,4,endpoint=False))
#fig5.set_size_inches(6, 3)

#### http://xkcd.com/color/rgb/

#sns.distplot( inner,fit=norm, kde = False, color = sns.xkcd_rgb["dark red"]   )
sns.distplot( inner,bins=200, kde = False, color = "red",label = r'$\mathbf{\alpha_{intra}}(intra)$')
#sns.distplot( inner,  kde = False, color = "red"   )

(mu, sigma) = norm.fit(inner)
print "mu={0}, sigma={1}".format(mu,sigma)

#sns.distplot( inter,fit=norm, kde = False, color = sns.xkcd_rgb["marine"]   )
sns.distplot( inter,bins=200, kde = False, color = "blue", label = r'$\mathbf{\alpha_{inter}}(inter)$'  )
#sns.distplot( inter,  kde = False,color = "blue"   )

(mu, sigma) = norm.fit(inter)
print "mu={0}, sigma={1}".format(mu,sigma)

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

plt.xlabel('Rotaiton Angle ($^\circ$)')
plt.ylabel('Probability Density')

plt.legend(loc=2)


#plt.legend(loc='best')
fig5.savefig('inner-inter-distribution-no-border.pdf',bbox_inches='tight')

#
#mu=20.6233143833, sigma=12.0882317208
#mu=-23.4042820655, sigma=13.8784374039

"""
