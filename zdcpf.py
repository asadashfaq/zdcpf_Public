#! /usr/bin/env python
from pylab import *
from scipy import *
import scipy.optimize as optimize
from scipy.sparse import coo_matrix
from numpy import concatenate as conc
import numpy as np
from time import time
import sys, os
from copy import deepcopy
import ctypes as ct

colors_countries = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers.

def get_positive(x):
    return x*(x>0.)  #Possibly it has to be x>1e-10.

class ESS:
    def __init__(self,pmax,emax,soco,nhours,inflow=0):
        self.SOC=np.zeros(nhours)
        self.power=np.zeros(nhours)
        self.SOC[0]=soco
        self.pmax=pmax
        self.pmin=-pmax
        self.emax=emax

    def avail(self,t):
        pout = min(self.pmax,self.SOC[t]*self.emax)       ## least of Discharge capacity and available Energy
        pin = min(-self.pmin,(1-self.SOC[t])*self.emax)   ## least of Charge capacity and available space
        return pout,pin

    def change(self,delta,t):
        #if delta > 0:  ## Storage will provide power, we are discharging
        self.SOC[t+1] = self.SOC[t] - delta/self.emax
        if self.SOC[t+1]>1: self.SOC[t+1]=1
        if self.SOC[t+1]<0: self.SOC[t+1]=0
        # if delta < 0:  ## Storage will absorb, we are charging
        ## but then command is the same... SOC -= delta/self.emax

class node:
    def __init__(self,path,fileName,ID):
        self.id = ID
        data = load(path + fileName)
        self.gamma = 1.0#float(setup[ID][0])
        self.alpha = 0.725#float(setup[ID][1]) #Alpha should be expanded to a vector. completalpha() can be applied in update()
        self.load = 1000*array(map(double,data['L']))
        self.nhours = len(self.load)
        self.normwind = array(map(double,data['Gw']))
        self.normsolar = array(map(double,data['Gs']))
        self.mean = mean(self.load)
        self.balancing = np.zeros(self.nhours)
        self.curtailment = np.zeros(self.nhours)
        self.label = data['datalabel']
        self.mismatch = None
        self.colored_import = None #Set using self.set_colored_i_import()
        data.close()
        self._update_()
        self.gen=np.zeros(self.nhours)


    def _update_(self):
        self.mismatch=(self.get_wind()+self.get_solar())-self.load
    
    def get_import(self):
        """Returns import power time series in units of MW."""
        return get_positive(get_positive(-self.mismatch) - self.balancing) #Balancing is exported if it exceeds the local residual load.
        
    def get_export(self):
        """Returns export power time series in units of MW."""
        return get_positive(self.mismatch) - self.curtailment #+ get_positive(self.balancing - get_positive(-self.mismatch))

    def get_localRES(self):
        """Returns the local use of RES power time series in units of MW."""
        return self.get_wind() + self.get_solar() - self.curtailment - self.get_export()

    def get_localBalancing(self):
        """Returns the local use of balancing power time series in units of MW."""
        return get_positive(-self.mismatch) - self.get_import()

    def get_wind(self):
        """Returns wind power time series in units of MW."""
        return self.mean*self.gamma*self.alpha*self.normwind

    def get_solar(self):
        """Returns solar power time series in units of MW."""
        return self.mean*self.gamma*(1.-self.alpha)*self.normsolar

    def set_gamma(self,gamma,operation='='):
        if operation == '=':
            self.gamma = gamma
        else:
            self.gamma *= gamma
        self._update_()
    
    def set_alpha(self,alpha):
        self.alpha=alpha
        self._update_()



class Nodes:
    def __init__(self,path='./data/',files=['ISET_country_AT.npz','ISET_country_FI.npz','ISET_country_NL.npz','ISET_country_BA.npz','ISET_country_FR.npz','ISET_country_NO.npz','ISET_country_BE.npz','ISET_country_GB.npz','ISET_country_PL.npz','ISET_country_BG.npz','ISET_country_GR.npz','ISET_country_PT.npz','ISET_country_CH.npz','ISET_country_HR.npz','ISET_country_RO.npz','ISET_country_CZ.npz','ISET_country_HU.npz','ISET_country_RS.npz','ISET_country_DE.npz','ISET_country_IE.npz','ISET_country_SE.npz','ISET_country_DK.npz','ISET_country_IT.npz','ISET_country_SI.npz','ISET_country_ES.npz','ISET_country_LU.npz','ISET_country_SK.npz'],load_filename=None):
        self.cache=[]
        for i in range(len(files)):
            n=node(path,files[i],i)
            self.cache=append(self.cache,n)
        F=np.zeros((size(files),self.cache[0].nhours))
        
        if load_filename != None:
            self._load_nodes_(load_filename,path='./results/')

    def __getitem__(self,x):
        return self.cache[x]
        
    def __len__(self):
        return len(self.cache)

    def set_gammas(self,value):
        # to change a single node's gamma, just write XX.set_gamma(yy)
        # 'value' can be a single number or a vector
        if np.size(value)==1:
            for i in self.cache: i.set_gamma(value)
        elif np.size(value)!=np.size(self.cache):
            print "Wrong gamma vector size. ", np.size(value,0)," were received, ",np.size(self.cache)," were expected."
        else:
            for i in self.cache:
                i.set_gamma(value[i.id])

    def set_alphas(self,value):
        # to change a single node's alpha, just write XX.set_gamma(yy)
        # 'value' can be a single number or a vector
        if size(value)==1:
            for i in self.cache: i.set_alpha(value)
        elif size(value)!=size(self.cache):
            print "Wrong gamma vector size. ", size(value,0)," were received, ",size(self.cache)," were expected."
        else:
            for i in self.cache:
                i.set_alpha(value[i.id])

    def save_nodes(self,filename,path='./results/'):
        """Saves the contents of a Nodes instance to a npz file."""
        
        attribute = dir(self[0])
        save_str = []
        #Determine which attributes to be saved
        for attribute in dir(self[0]):
            if attribute[0]=='_':
                continue
            elif is_numlike(getattr(self[0],attribute)) or is_string_like(getattr(self[0],attribute)):
                save_str.append(attribute + '=' + 'array([self[i].'+attribute+' for i in arange(len(self))])')

        #Write save file
        eval('np.savez(path+filename,'+','.join(save_str)+')')

        print 'Saved nodes to file: ', path+filename
        sys.stdout.flush()
        
    def _load_nodes_(self,load_filename,path='./results/'):
        """Loads a Nodes instance from an npz file."""

        npzobj = np.load(path+load_filename)
        
        for attribute in npzobj.files: ## this is real
            if (attribute == 'balancing') or (attribute =='curtailment'): ## this is my bullshit
                for i in arange(len(self)):
                    setattr(self.cache[i],attribute,npzobj[attribute][i])

# for i in arange(len(self)):
# print self.cache[i],attribute,npzobj[attribute][i]
# setattr(self.cache[i],attribute,npzobj[attribute][i])
        
        for n in self.cache:
            n._update_()
            
        print 'Loaded nodes from file: ', path+load_filename
        sys.stdout.flush()



def AtoKh(N,pathadmat='./settings/admat.txt'):
    Ad=np.genfromtxt(pathadmat,dtype='d')
    L=0
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    L+=1
    K_values=[]
    K_column_indices=[]
    K_row_indices=[]
    h=np.zeros(L*2)
    h=np.append(h,np.zeros(3*len(Ad)))
    L=0
    j=0
    i=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    K_values.extend([1,-1])
                    K_column_indices.extend([L,L])
                    K_row_indices.extend([j,i])
                    h[2*L]=Ad[i,j]
                    h[2*L+1]=Ad[j,i]
                    listFlows.append([str(N[j].label)+" to " +str(N[i].label), L])
                    L+=1
    #K=spmatrix(K_values,K_row_indices,K_column_indices)
    K=coo_matrix((K_values,(K_row_indices,K_column_indices)))    
    return K,K_values,h, listFlows

def get_quant(quant=0.99,filename='results/copper_flows.npy'):
    f=np.load(filename)
    flows=[]
    for i in f:
        flows.append(i*(i>0.))
        flows.append(-i*(i<0.))
    a=np.zeros(len(flows))
    b=np.zeros(len(flows))
    hs=np.zeros(len(flows))
    for i in range(len(flows)):
        a=hist(flows[i],cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i]=a[1][j]
                break
    return hs


def sdcpfold(admat='admat.txt',path='./settings/',copper=0,lapse=None,b=1,h0=None):
    N=Nodes()
    firststep=ct.CDLL('./balmin/libbalmin.dylib')
    secondstep=ct.CDLL('./flowmin/libflowmin.dylib')
    firststep.balmin.restype=ct.c_double
    secondstep.flowmin.restype=ct.c_int
    if copper == 1:
        firststep=ct.CDLL('./ubalmin/libubalmin.dylib')
        secondstep=ct.CDLL('./uflowmin/libuflowmin.dylib') 
        firststep.ubalmin.restype=ct.c_double
        secondstep.uflowmin.restype=ct.c_int    

    if lapse == None:
        lapse=N[0].nhours
    kv, Km, H,Lf=AtoKh(N)
    if h0 <> None:
        H=h0
    h_neg=-b*H[1:88:2]
    h_pos=b*H[0:88:2]
    k=array([float(i) for i in kv])
    Nlin=44
    Nnod=27  
   
    K=k.ctypes.data_as(ct.c_void_p)
    H_neg=h_neg.ctypes.data_as(ct.c_void_p)
    H_pos=h_pos.ctypes.data_as(ct.c_void_p)
    x=np.zeros(Nlin)
    flw=x.ctypes.data_as(ct.POINTER(ct.c_double))

    Flows=np.zeros((Nlin,lapse))
    start=time()
    delta=np.zeros(Nnod)
    g=np.zeros(Nlin)
    for t in range(lapse):
        for i in N:
            delta[i.id]=i.mismatch[t]
        Delta=delta.ctypes.data_as(ct.c_void_p)
        if copper == 0: MinBal=firststep.balmin(Delta,K,H_neg,H_pos)#.data_as(ct.c_double)
        if copper == 1: MinBal=firststep.ubalmin(Delta,K)#.data_as(ct.c_double)
        minbal=ct.c_double(MinBal+10.0)
        if copper == 0: MinFlow=secondstep.flowmin(Delta,K,H_neg,H_pos,minbal,flw)
        if copper == 1: MinFlow=secondstep.uflowmin(Delta,K,minbal,flw)
        #print MinBal
        for j in range(Nlin):
            Flows[j][t]=flw[j]

        if np.mod(t,1000)==0:
            print "t = ",t,". Elapsed time is ", time()-start
    print "Time series complete, elapsed time is ", round(time()-start)
    print "Attempting Nodal Assignment"
    Km=np.matrix(Km.todense())   # Turning the sparse matrix into a dense matrix
    for t in range(lapse):
        for i in N:
            delta[i.id] = i.mismatch[t]
        T=delta-(Km*Flows[:,t])#dot(Km,Flows[:,t])  # unrefined bal/curt  ===   Delta - KF
        balancing = -1*T*(1-sign(T))/2.0
        curtailment = T*(1+sign(T))/2.0
        for i in N:
            i.balancing[t] = balancing[i.id]
            i.curtailment[t] = curtailment[i.id]
    end=time()
    print "Nodes updated, Elapsed time is", round(end-start)
    return N,Flows       

def sdcpf(N,admat='admat.txt',path='./settings/',copper=0,lapse=None,b=1.0,h0=None):
    Nlinks=44
    Nnodes=27
    firststep=ct.CDLL('./balmin/libbalmin.dylib')
    firststep.balmin.restype=ct.c_double # default return type is int
    ufirststep=ct.CDLL('./ubalmin/libubalmin.dylib')
    ufirststep.ubalmin.restype=ct.c_double # default return type is int
    secondstep=ct.CDLL('./flowmin/libflowmin.dylib')
    secondstep.flowmin.restype=ct.c_int # just for fun
    usecondstep=ct.CDLL('./uflowmin/libuflowmin.dylib')
    usecondstep.uflowmin.restype=ct.c_int # just for fun

    if lapse == None:
        lapse=N[0].nhours
    km,kv,H,Lf=AtoKh(N) # dummy node has been deleted from admat.txt!!!
    km = np.matrix(km.todense())
    if (h0 != None):
        H=h0
    h_neg=b*-H[1:88:2]
    h_pos=b*H[0:88:2]
    #print 'h_pos: ',shape(h_neg),h_neg
    if (copper == 1):
        h_neg=-1.e6*np.ones(Nlinks)
        h_pos=1.e6*np.ones(Nlinks)
    #print 'h_pos: ',shape(h_neg),h_neg
    flw=np.zeros(Nlinks)
    k=array([float(i) for i in kv])
    K=k.ctypes.data_as(ct.c_void_p)
    H_neg=h_neg.ctypes.data_as(ct.c_void_p)
    H_pos=h_pos.ctypes.data_as(ct.c_void_p)
    Flw=flw.ctypes.data_as(ct.POINTER(ct.c_double)) # for return by ref

    delta=np.zeros(Nnodes)

    F=np.zeros((Nlinks,lapse)) # save flows and deltas to calc bal and curt later
    deltas=np.zeros((Nnodes,lapse))
    eps=1e-1
    start=time()
    if (copper == 0):
        for t in range(lapse):
            for i in N:
                delta[i.id]=i.mismatch[t]
            deltas[:,t]=delta
            Delta=delta.ctypes.data_as(ct.c_void_p)
            MinBal=firststep.balmin(Delta,K,H_neg,H_pos)
            # print "MinBal is ", MinBal
            minbal=ct.c_double(MinBal+eps)
            dummy=secondstep.flowmin(Delta,K,H_neg,H_pos,minbal,Flw)
            for i in range(Nlinks):
                F[i,t]=Flw[i]
            # print "MinFlows are "
            # print F[:,t]
            end=time()
            if (np.mod(t,2073)==0) and t>0:
                 print "Elapsed time is %3.1f seconds. t = %u out of %u" % ((end-start), t, lapse)
                 sys.stdout.flush()
        end=time()
        print "Calculation took %3.1f seconds." % (end-start)
        sys.stdout.flush()
    else: # use special unbounded C-fctns for unbounded copper flow
        for t in range(lapse):
            for i in N:
                delta[i.id]=i.mismatch[t]
            deltas[:,t]=delta
            Delta=delta.ctypes.data_as(ct.c_void_p)
            MinBal=ufirststep.ubalmin(Delta,K)
            # print "MinBal is ", MinBal
            minbal=ct.c_double(MinBal+eps)
            dummy=usecondstep.uflowmin(Delta,K,minbal,Flw)
            for i in range(Nlinks):
                F[i,t]=Flw[i]
            # print "MinFlows are "
            # print F[:,t]
            end=time()
            if (np.mod(t,2073)==0) and t>0:
                print "Elapsed time is %3.1f seconds. t = %u out of %u" % ((end-start), t, lapse)
                sys.stdout.flush()
        end=time()
        print "Calculation took %3.1f seconds." % (end-start)
        sys.stdout.flush()


    # no matter how the flows were obtained, we still have to calc bal and curt
    start2=time()
    for t in range(lapse):
        tmp=np.array(F[:,t])
        tmp2=np.array(deltas[:,t])
        tmp3=transpose(tmp2-dot(km,tmp))
        balancing=[-tmp3[i,0] if tmp3[i,0]<0 else 0. for i in range(len(tmp3))]
        curtailment=[tmp3[i,0] if tmp3[i,0]>0 else 0. for i in range(len(tmp3))]
        for i in N:
            i.balancing[t] = balancing[i.id]
            i.curtailment[t] = curtailment[i.id]
    end=time()
    print "Assigning balancing and curtailment took %3.1f seconds." % (end-start2)
    print "Complete calculation took %3.1f seconds." % (end-start)
    return N,F


#End of Code-Code

##




