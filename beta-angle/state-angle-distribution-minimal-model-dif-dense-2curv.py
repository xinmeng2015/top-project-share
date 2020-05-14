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


from sklearn import mixture
import matplotlib.pyplot
import matplotlib.mlab




InFile = 'inner-inter-state-dense.txt' #'states-angle-2d-noborder.txt'


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
 
for j in xrange(196000):
    demo = mix[j*2] - mix[j*2+1]
    dif.append(demo)


### plot distribution


# the histogram of the data
n, bins, patches = plt.hist(dif, 100, normed=1)

xx = []

for i in xrange(100):
    xx.append([bins[i],n[i]])

#print xx

#data = np.array(xx)

demo = np.array(dif)
data = demo[:,np.newaxis]

#print data


#print (n)
#print (bins)
#print (patches)


fig5, ax5 = plt.subplots()




###

clf = mixture.GMM(n_components=2, covariance_type='full')


#data = np.array(dif)
clf.fit(data)


m1, m2 = clf.means_
w1, w2 = clf.weights_
c1, c2 = clf.covars_

print 'm1',m1,m2
print 'w1',w1,w2
print 'c1',c1,c2


histdist = plt.hist(dif, 100, normed=True)
plotgauss1 = lambda x: plt.plot(x,w1*matplotlib.mlab.normpdf(x,m1,np.sqrt(c1))[0], linewidth=3)
plotgauss2 = lambda x: plt.plot(x,w2*matplotlib.mlab.normpdf(x,m2,np.sqrt(c2))[0], linewidth=3)
plotgauss1(histdist[1])
plotgauss2(histdist[1])

plotgauss3 = lambda x: plt.plot(x,w2*matplotlib.mlab.normpdf(x,m2,np.sqrt(c2))[0] + w1*matplotlib.mlab.normpdf(x,m1,np.sqrt(c1))[0] , linewidth=3)
plotgauss3(histdist[1])

#ax5.set_xlim(-70.0,70.0)
#ax5.set_xticks(np.linspace(-60,60,7,endpoint=True))
#ax5.set_ylim(0,0.04)
#ax5.set_yticks(np.linspace(0,0.04,4,endpoint=False))
#fig5.set_size_inches(6, 3)

#### http://xkcd.com/color/rgb/

#sns.distplot( dif, kde = True, color = "red",label = 'inter-intra-minimal-model'   )
#sns.distplot( dif,bins=100, kde = False, color = "red",label = 'inter-intra-minimal-model'   )
#sns.distplot( dif, fit=norm, kde = False, color = "red",label = 'inter-intra-minimal-model'   )

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

plt.xlabel('Rotaiton Angle ($^\circ$)')
plt.ylabel('Probability Density')

plt.legend(loc='best')


#plt.legend(loc='best')
fig5.savefig('minimal-model-inner-inter-dif-dense-fit.pdf',bbox_inches='tight')

#
#mu=19.9064230234, sigma=14.5756637612
#mu=-22.9004333214, sigma=16.4392568036


