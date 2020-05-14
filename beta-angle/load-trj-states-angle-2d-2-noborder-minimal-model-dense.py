#!/usr/bin/:wenv python
#_*_ coding: utf-8 _*_

__doc__ = """
/****************************************************************************************************************************
                                  SQN-PYTHON: call gromacs to calculate force
Catogory:   testing
Name:   1-test.py
Previous Version:   6-test.py
Description && Note:
                    read gromacs trr file
                    and do coordinates based analysis

                    loadGro()  last z rows should be 36:44
                    now can calculate the whole strcuture, do not search the whole
History:
        <Start Date>:   20-07-2016
        <Editing Date>: *-07-2016
        <Finish Date>:
*********************************************************************************************/
"""

#import pyopencl as cl
import numpy as np
import math  as m
from scipy import linalg
import matplotlib.pyplot as plt
import sh
#import MDAnalysis as mda
import MDAnalysis.coordinates.TRR as trr
import MDAnalysis.coordinates.XTC as xtc
#import MDAnalysis.coordinates.TRR.TRRReader as TRRreader
####
###
InTrrFile  = 'traj_comp.xtc' #'traj.trr'
### 1 frame information based on topol.top
### topol.top: MOL 50
### mol.itp: atoms 121 , Mg_index=27(start from 1)
N_PAIR      = 900 ## 28 * 28
N_MOL       = 1800 ## 28 * 28 * 2
MOL_NATOM   = 121
Mg_index    = 27 ## start from 1
Carbonyl_O_index  = 56 #start from 1
Hydroxyl_O_index  = 33 #start from 1
Hydroxyl_H_index  = 34 #start from 1
N_R1_index  = 3 #start from 1
N_R2_index  = 11 #start from 1
N_R3_index  = 17 #start from 1

####
MG_OH_coupling = 0.3 ## criteria for couping between MG and OH
rDA = 0.35 ## hydrogen bond geometry criteria
alphaHDA = 30

###################################################
packFile =  'inner-inter-state-dense.txt' #'states-angle-2d-no-border-minimal-model.txt'
PackAFile= 'packA.txt'
PackBFile= 'packB.txt'
packGammaFile= 'packGamma.txt'

####
InGroFile = 'derote-confout.gro'
REMARK = 2
ATOMNUM = 217800 #11858

SEARCH_CONSTAINTS = 80


# A --- O2
# O --- Mg
# B --- CR1R2     carbon connect ring1 and ring2
# !! calculate crossover C AND D
# calculate BOA  and INNER INTER  

A_index = 56  ## O2
O_index = 27  ## Mg
B_index = 6    # CR1R2


#####
def calc_dist(p1,p2):
    return m.sqrt((p1[0] - p2[0]) ** 2 +(p1[1] - p2[1]) ** 2 +(p1[2] - p2[2]) ** 2)
    #return m.sqrt((p2[0] - p1[0]) ** 2 +(p2[1] - p1[1]) ** 2 +(p2[2] - p1[2]) ** 2)
def calc_dist_1(p1,p2):
    return m.sqrt((m.fabs(p1[0] - p2[0])+CLAW) ** 2 +(p1[1] - p2[1]) ** 2 +(p1[2] - p2[2]) ** 2)
def calc_dist_2(p1,p2):
    return m.sqrt((p1[0] - p2[0]) ** 2 +(m.fabs(p1[1] - p2[1])+CLAW) ** 2 +(p1[2] - p2[2]) ** 2)
def calc_dist_3(p1,p2):
    return m.sqrt((p1[0] - p2[0]) ** 2 +(p1[1] - p2[1]) ** 2 + (m.fabs(p1[2] - p2[2])+CLAW) ** 2)
#################################################################################################


def call_gromacs(network):

    trj = xtc.XTCReader(InTrrFile) ## trr.TRRReader


    ### test reading first and last frame and compare with gro file is ok
    """
    for t in xrange(1): ## first frame test
        for i in xrange(N_MOL):
            #if ((i+1)% 2 == 1): ## the
            #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0 ## unit A-->nm
                #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ N_R1_index-1]/10.0
                #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ N_R2_index-1]/10.0
                #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ N_R3_index-1]/10.0
            print "gmx trr position",trj[1000].positions[i*MOL_NATOM+ Mg_index-1]/10.0 ## unit A-->nm
    """
    mgmg_dist_list = []
    state_list = []

    l0_all = []
    l1_all = []
    aob_all = []
    inner_raw_all = []
    inter_raw_all = []
    inner_norm_all = []
    inter_norm_all = []
    l2_all = []
    l3_all = []
    l2_2_all = []
    l3_2_all = []




    for t in xrange(1000): ## first frame test
        #print t
        if (t > 500 and (t+1) % 2 ==0):
            ##inside frame t the angle calculation, just need the
            for i in xrange(60,N_MOL-60):
                if (i % 60 > 1 and i % 60 < 58 ):
                    print t, i
                    if (network[i][0] != -1):
                        target    =  int(network[i][0])

                        demoo1 = i*MOL_NATOM+ O_index  -1
                        o1 = np.copy ( trj[t].positions[demoo1]/10.0 )
                        demoo2 = target*MOL_NATOM+ O_index  -1
                        o2 = np.copy ( trj[t].positions[demoo2]/10.0 )

                        #demoa1 = i*MOL_NATOM+ A_index  -1
                        #a1 = np.copy ( trj[t].positions[demoa1]/10.0 )
                        demoa2 = target*MOL_NATOM+ A_index  -1
                        a2 = np.copy ( trj[t].positions[demoa2]/10.0 )

                        demob1 = i*MOL_NATOM+ B_index  -1
                        b1 = np.copy ( trj[t].positions[demob1]/10.0 )
                        #demob2 = target*MOL_NATOM+ B_index  -1
                        #b2 = np.copy ( trj[t].positions[demob2]/10.0 ) 


                        v_ob1 = ( b1 - o1 )
                        v_oa2 = ( a2 - o2 )

                    if (i%2 == 0):
                        cosinner = np.dot(v_ob1,v_oa2)/(linalg.norm(v_ob1)*linalg.norm(v_oa2))
                        inner = np.rad2deg ( np.arccos( cosinner ))
                        inner_raw_all.append(inner )
                        state_list.append(State(i+1,target+1,inner))
                    else:
                        cosinter = np.dot(v_ob1,v_oa2)/(linalg.norm(v_ob1)*linalg.norm(v_oa2))
                        inter = np.rad2deg ( np.arccos( cosinter ))
                        inter_raw_all.append(inter )
                        state_list.append(State(i+1,target+1,inter))

                        """
                    	Mgi = np.copy ( trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0 )
                    	Mgj = np.copy ( trj[t].positions[target*MOL_NATOM+ Mg_index-1]/10.0 )
                    	dist0 = calc_dist(Mgj, Mgi)
                    	mgmg_dist_list.append( dist0 )

                    	demoi1 = np.copy( trj[t].positions[i*MOL_NATOM+ N_R1_index-1]/10.0 )
                    	demoi2 = np.copy( trj[t].positions[i*MOL_NATOM+ N_R2_index-1]/10.0 )
                    	demoi3 = np.copy( trj[t].positions[i*MOL_NATOM+ N_R3_index-1]/10.0 )
                    	vi1 = (demoi1-demoi2)
                    	vi2 = (demoi3-demoi2)
                    	
                    	## 1 i surface norm vector
                    	normalVectori = np.cross(vi1,vi2)/linalg.norm(np.cross(vi1,vi2))
                    	## 2 rotate part.
                    	## translation to make a origin point
                    	## here we select the i Mg atom as the O point.
                    	## select origin point
                    	d = np.copy(trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0)
                    	Mgi -= d
                    	Mgj -= d
                    	O2i = np.copy(trj[t].positions[i*MOL_NATOM+ Carbonyl_O_index-1]/10.0 )
                    	O2i -= d
                    	
                    	a = normalVectori  ##normalized vector a
                    	
                    	#print "!!!!",a, linalg.norm(a)
                    	#print a[2]
                    	if a[2] > 0:
                    	    b = np.array([0.0,0.0,1.0])
                    	else:
                    	    b = np.array([0.0,0.0,-1.0]) ## make it along z direction
                    	### when do the rotation, syn and anti rotate to a same configuration so there is no different now!!!!
                    	### the judge is really important
                    	### show what the order one packing angle or the minimal model delta angle
                    	
                    	v = np.cross(a,b)     # a X b
                    	s = linalg.norm(v)
                    	c = np.dot(a,b)
                    	
                    	v_skew = np.array([[0.0,-v[2], v[1]],[v[2],0.0,-v[0]],[-v[1],v[0],0.0 ]])
                    	R = np.identity(3) + v_skew + np.dot(v_skew,v_skew)* (1.0-c) / (s*s)
                    	
                    	#print "mgi mgj  dist", Mgi, Mgj, calc_dist(Mgi,Mgj)
                    	
                    	demo_v1 = Mgj - Mgi
                    	demo_v2 = O2i - Mgi
                    	
                    	vectorMgMg = np.dot(R,demo_v1)
                    	vectorMgO2 = np.dot(R,demo_v2)
                    	
                    	vectorx = np.array( [    vectorMgO2[0], vectorMgO2[1] ]   )
                    	vectory = np.array( [    vectorMgMg[0], vectorMgMg[1] ]   )
                    	
                    	
                    	cosx = np.dot(vectorx ,vectory)/(linalg.norm(vectorx)*linalg.norm(vectory))
                    	angle = np.rad2deg ( np.arccos( cosx  ))
                    	
                    	sinx = np.cross(vectorx ,vectory)/(linalg.norm(vectorx)*linalg.norm(vectory))
                    	
                    	if sinx > 0:
                    	    sign = 1
                    	else:
                    	    sign = -1
                    	
                    	angle22 = angle * sign
                    	#print "--->", angle2, angle
                    	
                    	state_list.append( State(i+1,target+1,angle22))
                        """
    with open(packFile,'w') as fo:
        for x in state_list:
             fo.write("%d   %d   %f  "% ( x.indexi,x.indexj,x.angle ) + '\n')
             ## difference index is 1 inter dimer > 1 is inter dimer
             
"""
            Mg_dist=[]
            MgOH_dist=[]
            OHMg_dist=[]
            DA_dist =[]
            HDA_angle = []
            judge_Mg_OH = np.zeros(N_MOL)
            judge_OH_Mg = np.zeros(N_MOL)
            judge_HO_O2 = np.zeros(N_MOL)
            judge_O2_HO = np.zeros(N_MOL)
            judge_hbond = np.zeros(N_MOL)
            network     = np.zeros((N_MOL,2))-1
            #print network

            MgDistance_Matrix = np.zeros((N_MOL,N_MOL))

            packA     = []
            packB     = []
            packGamma = []

            for i in xrange(N_MOL-1):
                for j in xrange(i+1,N_MOL):

                    #print t,i,j  ## see loop

                    ## calculate Mg-Mg distances ## no periodic consideration
                    #demo1 = trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0
                    #demo2 = trj[t].positions[j*MOL_NATOM+ Mg_index-1]/10.0
                    #demo_dist = calc_dist(demo2, demo1)
                    #print demo_dist

                    ## calculate the stacking angle of two pigment ring
                    #demoi1 = trj[t].positions[i*MOL_NATOM+ N_R1_index-1]/10.0
                    #demoi2 = trj[t].positions[i*MOL_NATOM+ N_R2_index-1]/10.0
                    #demoi3 = trj[t].positions[i*MOL_NATOM+ N_R3_index-1]/10.0

                    #demoj1 = trj[t].positions[j*MOL_NATOM+ N_R1_index-1]/10.0
                    #demoj2 = trj[t].positions[j*MOL_NATOM+ N_R2_index-1]/10.0
                    #demoj3 = trj[t].positions[j*MOL_NATOM+ N_R3_index-1]/10.0

                    #print demoj3, demoj3.dtype
                    #vi1 = (demoi1-demoi2)
                    #vi2 = (demoi3-demoi2)
                    #print vi1, vi2
                    #vj1 = (demoj1-demoj2)
                    #vj2 = (demoj3-demoj2)                #print vi1,vj1,vi2,vj2
                    #normalVectori = np.cross(vi1,vi2)/linalg.norm(np.cross(vi1,vi2))
                    #normalVectorj = np.cross(vj1,vj2)/linalg.norm(np.cross(vj1,vj2))
                    #print normalVectori, normalVectorj, linalg.norm(normalVectori),linalg.norm(normalVectorj)
                    #cosij=np.dot(normalVectori,normalVectorj)/(linalg.norm(normalVectori)*linalg.norm(normalVectorj))
                    #print i,j,cosij
                    ## defect flip of the N1,N2,N3, Sin-Anti difference
                    ######### not useful
                    #angle0 = np.rad2deg (np.arccos(cosij))
                    #angle0 = np.rad2deg ( np.arccos( np.dot(normalVectori,normalVectorj) ) )
                    #print i,j,angle0  ## this is the stacking angle of two planes
                    #vectorN2jN2i = (demoj2-demoi2)
                    #vectorMgMg = (demo2-demo1)
                    #cosMgNorm = np.dot(normalVectori,vectorMgMg)/(linalg.norm(normalVectori)*linalg.norm(vectorMgMg))
                    #angle1 = np.rad2deg ( np.arccos( cosMgNorm  ))
                    #print angle1
                    #vectorN2N2 = (demoj2-demoi2)
                    #cosN2N2 = np.dot(normalVectori,vectorN2N2)/(linalg.norm(normalVectori)*linalg.norm(vectorN2N2))
                    #angle2 = np.rad2deg ( np.arccos( cosN2N2  ))
                    #print i,j,angle1, angle2
                    #############################
                #if ((i+1)% 2 == 1): ## the
                    #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0 ## unit A-->nm
                    #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ N_R1_index-1]/10.0
                    #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ N_R2_index-1]/10.0
                    #print "gmx trr position",trj[t].positions[i*MOL_NATOM+ N_R3_index-1]/10.0

                ## detect paired Mg distance
                    Mgi = trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0
                    Mgj = trj[t].positions[j*MOL_NATOM+ Mg_index-1]/10.0
                    dist0 = calc_dist(Mgj, Mgi)
                    # print i,j,dist0 correct
                    #Mg_dist.append(dist0)
                    MgDistance_Matrix[i][j]= dist0
                    MgDistance_Matrix[j][i]= dist0

                ## detect Mg-O(H) and (H)O-Mg distance judge for pair coupling
                    Mgi = trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0
                    OHj = trj[t].positions[j*MOL_NATOM+ Hydroxyl_O_index-1]/10.0
                    dist1 = calc_dist(Mgi, OHj)
                    #print i,j, dist1
                    if (dist1 < MG_OH_coupling):
                        MgOH_dist.append(dist1)
                        judge_Mg_OH[i]=1
                        judge_OH_Mg[j]=1
                        network[i][0] = j

                    Mgj = trj[t].positions[j*MOL_NATOM+ Mg_index-1]/10.0
                    OHi = trj[t].positions[i*MOL_NATOM+ Hydroxyl_O_index-1]/10.0
                    dist2 = calc_dist(Mgj, OHi)
                    #print i,j, dist2
                    if (dist2 < MG_OH_coupling):
                        OHMg_dist.append(dist2)
                        judge_OH_Mg[i]=1
                        judge_Mg_OH[j]=1
                        network[j][0] = i

                ##detect OH-02 /  O2-OH distance  and H-D --- A
                    OHi = trj[t].positions[i*MOL_NATOM+ Hydroxyl_O_index-1]/10.0
                    O2j = trj[t].positions[j*MOL_NATOM+ Carbonyl_O_index-1]/10.0
                    dist3 = calc_dist(OHi, O2j)
                    if (dist3 <= rDA):
                        DA_dist.append(dist3)
                        judge_HO_O2[i]=1
                        judge_O2_HO[j]=1
                        network[i][1] = j

                    HOi = trj[t].positions[i*MOL_NATOM+ Hydroxyl_H_index-1]/10.0
                    vDH = HOi - OHi
                    vDA = O2j - OHi
                    alpha0 = np.rad2deg( np.arccos( np.dot(vDH,vDA)/(linalg.norm(vDH)* linalg.norm(vDA))) )
                    #print i,j, alpha0  ## initial two kind of angles ~ 18 <30 and ~56 > 30
                    if (dist3 <= rDA and alpha0 <= alphaHDA):
                        HDA_angle.append(alpha0)
                        judge_hbond[i] = 1
                        judge_hbond[j] = 1
                    ##--
                    O2i = trj[t].positions[i*MOL_NATOM+ Carbonyl_O_index-1]/10.0
                    OHj = trj[t].positions[j*MOL_NATOM+ Hydroxyl_O_index-1]/10.0
                    dist4 = calc_dist(O2i, OHj)
                    if (dist4 <= rDA):
                        DA_dist.append(dist4)
                        judge_O2_HO[i]=1
                        judge_HO_O2[j]=1
                        network[j][1] =i

                    HOj = trj[t].positions[j*MOL_NATOM+ Hydroxyl_H_index-1]/10.0
                    vDH = HOj - OHj
                    vDA = O2i - OHj
                    alpha1 = np.rad2deg( np.arccos( np.dot(vDH,vDA)/(linalg.norm(vDH)* linalg.norm(vDA))) )
                    #print i,j, alpha0
                    if (dist4 <= rDA and alpha1 <= alphaHDA):
                        HDA_angle.append(alpha1)
                        judge_hbond[i] = 1
                        judge_hbond[j] = 1
                    ## Hydrogen bonds geometrical criterion ## gromacs 5 mannul 8.12 hydrogen bond
                    ## r(DA) <= 0.35 (canbe bigger like 0.4)
                    ## alpha(H-D-A) <= 30 degree


            ## start a b and gamma parameter check
            for i in xrange(N_MOL-1):
                if (network[i][0] != -1):
                    stack    =  network[i][0]
                    print "!!!!",i,stack
                    if (network[stack][0] != -1):
                        target1  =  network[stack][0]
                        dist5    =  MgDistance_Matrix[i][target1]
                        packA.append(dist5)
                        vectorA  = trj[t].positions[target1*MOL_NATOM+ Mg_index-1]/10.0-trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0
                        if (network[stack][1] != -1):
                            target2    = network[stack][1]
                            dist6    =  MgDistance_Matrix[i][target2]
                            packB.append(dist6)
                            vectorB  = trj[t].positions[target2*MOL_NATOM+ Mg_index-1]/10.0-trj[t].positions[i*MOL_NATOM+ Mg_index-1]/10.0
                            angle1   = np.rad2deg( np.arccos( np.dot(vectorA,vectorB)/(linalg.norm(vectorA)* linalg.norm(vectorB))) )
                            packGamma.append(angle1)
            #print packA, packB,packGamma
            #print np.mean(packA), np.std(packA)
            #print np.mean(packB), np.std(packB)
            #print np.mean(packGamma), np.std(packGamma)
            with open(packFile,'a') as fo:
                fo.write("%d        %f    %f    %f     %f     %f      %f" % (t, np.mean(packA),np.std(packA), np.mean(packB),np.std(packB),np.mean(packGamma), np.std(packGamma))+'\n')



            ##print DA_dist, HDA_angle, judge_hbond
            ##print MgOH_dist,OHMg_dist,judge_Mg_OH,judge_OH_Mg  ## boarder molecule has defect







                #print "gmx trr forces", trj[0].forces

        #self.vdphi = -np.reshape(trj[0].forces, NDIM)/4.184 ## unit from kj/mol/A  --> kcal/mol/A
        #print "gmx trr vdphi",self.vdphi


"""


############

def loadPosition(inputfilename):
    demo=[]
    GroAtom_list = []
    ##    1ETH     CB    1   1.460   1.559   1.491  0.2285  0.2711 -0.7051
    ##    0         1    2    3       4       5       6
    with open(inputfilename,'r') as f:
        data = f.readlines()
        for line in data[REMARK:REMARK+ATOMNUM]:
            #print (line)
            list = line.split()
            ## the resid and residuename belong to list[0],split it
            resid         = int(line[:5])
            residuename   = line[5:10]
            atomname      = line[10:15]
            index         = int(line[15:20])
            x             = float(line[20:28])  ## double in c, float64 in numpy --> doesn't matter in here
            y             = float(line[28:36])
            z             = float(line[36:44])
            vx            = 0.0 #float(list[6])
            vy            = 0.0 #float(list[7])
            vz            = 0.0 #float(list[8])
            #print x
            demo.append(x)
            demo.append(y)
            demo.append(z)
            GroAtom_list.append(GroAtom(resid,residuename, atomname, index, x, y, z,vx,vy,vz))
            #print x,y,z
    #print GroAtom_list
    #print "original", demo
    #self.vx=np.asarray(demo,dtype=np.float32) ## this float32 in needed, since in gpu is single precisoin
    #print "load  vx", self.vx, self.vx.shape, self.vx[0]

    return GroAtom_list  ## maybe self.GroAtom does not need to return


def calcPack(atoms):


    Mg_dist=[]
    MgOH_dist=[]
    OHMg_dist=[]
    DA_dist =[]
    HDA_angle = []
    judge_Mg_OH = np.zeros(N_MOL)
    judge_OH_Mg = np.zeros(N_MOL)
    judge_HO_O2 = np.zeros(N_MOL)
    judge_O2_HO = np.zeros(N_MOL)
    judge_hbond = np.zeros(N_MOL)
    network     = np.zeros((N_MOL,2))-1
    #print network

    MgDistance_Matrix = np.zeros((N_MOL,N_MOL))

    packA     = []
    packB     = []
    packGamma = []

    demox = [] #np.zeros []  np.zeros((2, 3)) a
    for atom in atoms:
        demox.append(atom.x)
        demox.append(atom.y)
        demox.append(atom.z)
    trj = np.reshape(demox,(ATOMNUM,3))
    #print trj
    ########### need the resid index is correct continuous sequence?


    for i in xrange(N_MOL-1):
        print i
        for j in xrange(i+1,min(i+SEARCH_CONSTAINTS,N_MOL)):

        ## detect paired Mg distance
            Mgi = trj[i*MOL_NATOM+ Mg_index-1]
            Mgj = trj[j*MOL_NATOM+ Mg_index-1]
            dist0 = calc_dist(Mgj, Mgi)
            ### print i,j,dist0  ## correct
            #Mg_dist.append(dist0)
            MgDistance_Matrix[i][j]= dist0
            MgDistance_Matrix[j][i]= dist0

        ## detect Mg-O(H) and (H)O-Mg distance judge for pair coupling
            Mgi = trj[i*MOL_NATOM+ Mg_index-1]
            OHj = trj[j*MOL_NATOM+ Hydroxyl_O_index-1]
            dist1 = calc_dist(Mgi, OHj)
            #print i,j, dist1
            if (dist1 < MG_OH_coupling):
                MgOH_dist.append(dist1)
                judge_Mg_OH[i]=1
                judge_OH_Mg[j]=1
                network[i][0] = j

            Mgj = trj[j*MOL_NATOM+ Mg_index-1]
            OHi = trj[i*MOL_NATOM+ Hydroxyl_O_index-1]
            dist2 = calc_dist(Mgj, OHi)
            #print i,j, dist2
            if (dist2 < MG_OH_coupling):
                OHMg_dist.append(dist2)
                judge_OH_Mg[i]=1
                judge_Mg_OH[j]=1
                network[j][0] = i

        ##detect OH-02 /  O2-OH distance  and H-D --- A
            OHi = trj[i*MOL_NATOM+ Hydroxyl_O_index-1]
            O2j = trj[j*MOL_NATOM+ Carbonyl_O_index-1]
            dist3 = calc_dist(OHi, O2j)
            if (dist3 <= rDA):
                DA_dist.append(dist3)
                judge_HO_O2[i]=1
                judge_O2_HO[j]=1
                network[i][1] = j

            HOi = trj[i*MOL_NATOM+ Hydroxyl_H_index-1]
            vDH = HOi - OHi
            vDA = O2j - OHi
            alpha0 = np.rad2deg( np.arccos( np.dot(vDH,vDA)/(linalg.norm(vDH)* linalg.norm(vDA))) )
            #print i,j, alpha0  ## initial two kind of angles ~ 18 <30 and ~56 > 30
            if (dist3 <= rDA and alpha0 <= alphaHDA):
                HDA_angle.append(alpha0)
                judge_hbond[i] = 1
                judge_hbond[j] = 1
            ##--
            O2i = trj[i*MOL_NATOM+ Carbonyl_O_index-1]
            OHj = trj[j*MOL_NATOM+ Hydroxyl_O_index-1]
            dist4 = calc_dist(O2i, OHj)
            if (dist4 <= rDA):
                DA_dist.append(dist4)
                judge_O2_HO[i]=1
                judge_HO_O2[j]=1
                network[j][1] =i

            HOj = trj[j*MOL_NATOM+ Hydroxyl_H_index-1]
            vDH = HOj - OHj
            vDA = O2i - OHj
            alpha1 = np.rad2deg( np.arccos( np.dot(vDH,vDA)/(linalg.norm(vDH)* linalg.norm(vDA))) )
            #print i,j, alpha0
            if (dist4 <= rDA and alpha1 <= alphaHDA):
                HDA_angle.append(alpha1)
                judge_hbond[i] = 1
                judge_hbond[j] = 1
            ## Hydrogen bonds geometrical criterion ## gromacs 5 mannul 8.12 hydrogen bond
            ## r(DA) <= 0.35 (canbe bigger like 0.4)
            ## alpha(H-D-A) <= 30 degree

    ## start a b and gamma parameter check
    for i in xrange(N_MOL-1):
        if (network[i][0] != -1):
            stack    =  int(network[i][0])
            ####print "!!!!",i+1,stack+1
            if (network[stack][0] != -1):
                target1  =  int(network[stack][0])
                dist5    =  MgDistance_Matrix[i][target1]
                packA.append(dist5)
                vectorA  = trj[target1*MOL_NATOM+ Mg_index-1]-trj[i*MOL_NATOM+ Mg_index-1]
                if (network[stack][1] != -1):
                    target2    = int(network[stack][1])
                    dist6    =  MgDistance_Matrix[i][target2]
                    packB.append(dist6)
                    vectorB  = trj[target2*MOL_NATOM+ Mg_index-1]-trj[i*MOL_NATOM+ Mg_index-1]
                    angle1   = np.rad2deg( np.arccos( np.dot(vectorA,vectorB)/(linalg.norm(vectorA)* linalg.norm(vectorB))) )
                    packGamma.append(angle1)

    #print packA, packB,packGamma
    print np.mean(packA), np.std(packA)
    print np.mean(packB), np.std(packB)
    print np.mean(packGamma), np.std(packGamma)

    return network



class Atoms:
    def __init__(self,index,x,y,z,mass,sigma,epsilon,charge):
        self.index     = index
        self.x         = x
        self.y         = y
        self.z         = z
        self.mass      = mass
        self.sigma     = sigma
        self.epsilon   = epsilon
        self.charge    = charge



class GroAtom:
    def __init__(self,resid,residuename, atomname, index, x, y, z,vx,vy,vz):
        self.resid        = resid
        self.residuename  = residuename
        self.atomname     = atomname
        self.index        = index
        self.x         = x
        self.y         = y
        self.z         = z
        self.vx        = vx
        self.vy        = vy
        self.vz        = vz

class State:
    def __init__(self,indexi,indexj,angle):
        self.indexi     = indexi
        self.indexj     = indexj
        self.angle    = angle




if __name__=="__main__":


    groAtoms = loadPosition(InGroFile)
    couple_network = calcPack(groAtoms)
    call_gromacs(couple_network)



    ### Intial load
    #groAtoms = test.loadPosition("ethanol.gro")
    #itpAtoms,itpBonds,itpAngles,itpDihedrals) = test.loadItp(InItpFile)
    #oplsAtoms = test.loadItpFF(InItpFile2)

    #test.popCorn(groAtoms)
    #test.execute(groAtoms)

    ### visualizaiton
    #print system_potential2
    #test.visualization()
    ###
