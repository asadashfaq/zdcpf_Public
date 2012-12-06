#! /usr/bin/env python
from zdcpf import *
import networkx as nx
import csv

def AtoKh_nonsparse(N,pathadmat='./settings/admatold.txt'):
    Ad=np.genfromtxt(pathadmat,dtype='d')
    L=0
    G=nx.Graph()
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    L+=1
    K=np.zeros((len(Ad),L))
    h=np.zeros(L*2)
    h=np.append(h,np.zeros(3*len(Ad)))
    L=0
    j=0
    i=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    K[j,L]=1
                    K[i,L]=-1
                    h[2*L]=Ad[i,j]
                    h[2*L+1]=Ad[j,i]
                    if L>0: listFlows.append([str(N[j-1].label)+" to " +str(N[i-1].label), L-1])
                    G.add_edge(str(N[j-1].label),str(N[i-1].label),weight=max(h[2*L],h[2*L+1]))
                    L+=1
    return K,h,listFlows,G 

def show_hist(link,tit,filename='results/copper_flows.npy',e=1,b=500,limit=10000):
    f=np.load(filename)
    flows=[]
#    ax=subplot(1,1,1)
    for i in f:
        flows.append(i)
#        flows.append(i*(i>e))
#        flows.append(i*(i<-e))
#    a=hist(flows[link],bins=b,range=[0.1,flows[link].max()],normed=0,histtype='stepfilled')
#    b=hist(flows[link+1],bins=b,range=[flows[link+1].min(),-0.1],normed=0,histtype='stepfilled')
    a=hist(flows[link],bins=b,normed=0,histtype='stepfilled')
    plt.ylabel('Incidence')
    plt.xlabel('Transmission Magnitude (MW)')
    plt.title(r'Histogram of Power Flows : '+tit)
    vlines(-31602,0,150,color='r',linewidth=1.2)
    vlines(20265,0,150,color='r',linewidth=1.2)#,linestyle='dashed'
    plt.text(-37250,125,'1% Q')
    plt.text(20650,125,'99% Q')
#    plt.vlines(0,0,1000,linestyle='dashed')
#    plt.axvline(15357.73,0,800,linestyle='dashed',color='r')
#    plt.axvline(7678.86,0,800,linestyle='dashed',color='r')
#    plt.axvline(3993.00,0,800,linestyle='dashed',color='r')
#    plt.axvline(307.15,0,800,linestyle='dashed',color='r')
#    plt.text(15457,400,'99Q')
#    plt.text(7778,500,'95Q')
#    plt.text(4093,600,'90Q')
#    plt.text(407,700,'75Q')
    gcf().set_size_inches([9*1.5,3*1.75])
    plt.axis([-80000.1,40000.1,0,600])#plt.axis([0.5, 1, 0, 38000])
    plt.grid(True)
    savefig('./figures/' + str(link) +'.pdf', dpi=300)
    #show()


def savehist(link,name):
    f=np.load('results/copper_flows.npy')
    N=Nodes()
    a,b,c,d=AtoKh_nonsparse(N)
    a=hist(f[link],bins=500,normed=1)
    b=[a[1][:-1],a[0]]
    save('./figures/'+c[link][0],b)
#    g=csv.writer(open(name,'wb'),delimiter=',')
#    for column in f[link]: 
#        g.writerow(column)

def get_transmissions():
    NA=Nodes(load_filename='Case_A_Beta_1.0.npz')
    a,b,ha,d=AtoKh(NA)
    hc=get_quant(0.9999)
    hq=get_quant(0.99)
    hd=0.4*get_quant(0.99)

    HC=biggestpair(hc)
    HQ=biggestpair(hq)
    HD=biggestpair(hd)
    HA=biggestpair(ha)
    i=0
    while i<(44):
        print d[i][0],'&',str(round(hc[2*i]/1000,2)),'&',str(round(hq[2*i]/1000,2)),'&',str(round(hd[2*i]/1000,2)),'&',str(round(ha[2*i]/1000,2)),'\\\\'
        print '$\\blacktriangleleft$','&',str(round(hc[2*i+1]/1000,2)),'&',str(round(hq[2*i+1]/1000,2)),'&',str(round(hd[2*i+1]/1000,2)),'&',str(round(ha[2*i+1]/1000,2)),'\\\\'   
        i+=1 
    print 'Total & ', sum(HC),'&', sum(HQ),'&', sum(HD),'&', sum(HA),'\\\\'


def importexport():
    NC=Nodes(load_filename='copper_nodes.npz')
    NQ=Nodes(load_filename='Case_C_Beta_1.0.npz')
    ND=Nodes(load_filename='Case_C_Beta_0.4.npz')
    NA=Nodes(load_filename='Case_A_Beta_1.0.npz')
    impwanted=np.ones(28)
    imp=np.ones(28)
    i=[]
    imports=[]
    expwanted=np.ones(28)
    expo=np.ones(28)
    e=[]
    exports=[]
    CASES=[NC,NQ,ND,NA]
    imports=np.zeros((28,4))
    exports=np.zeros((28,4))
    for N in CASES:
        impwanted[27]=0
        expwanted[27]=0
        for n in N:
            impwanted[n.id]=sum(get_positive(-1*n.mismatch))
            expwanted[n.id]=sum(get_positive(n.mismatch))
            imp[n.id]=impwanted[n.id]-sum(n.balancing)
            #imp[27]=sum(imp[:-1])
            expo[n.id]=expwanted[n.id]-sum(n.curtailment)
            #expo[27]=sum(expo[:-1])
            x=imp[n.id]/impwanted[n.id]
            y=expo[n.id]/expwanted[n.id]
            #if n.id == 19:
            #    print impwanted,imp,x
            i.append(x)
            e.append(y)
        impwanted[27] = sum(impwanted[:-1])
        expwanted[27] = sum(expwanted[:-1])
        imp[27]=sum(imp[:-1])
        expo[27]=sum(imp[:-1])
        i.append(imp[27]/impwanted[27])
        e.append(expo[27]/expwanted[27])
        if N==NC: p=0    
        if N==NQ: p=1
        if N==ND: p=2
        if N==NA: p=3
        imports[:,p]=i[p*28:(p+1)*28]
        exports[:,p]=e[p*28:(p+1)*28]       

    return imports,exports


#scenario A
def Case_A(betas=[0.0,0.25,0.5,0.75,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0]):
    N=Nodes()    
    for b in betas:
        N,F=sdcpf(N,b=b)
        N.save_nodes('Case_A_Beta_'+str(b))
        save('./results/'+'Flows_Case_A_Beta_'+str(b),F)

def Case_B(links=np.arange(0.0,30000.1,1500.0)):
    hopt=get_quant(.99)
    h0=get_quant(.99)
    N=Nodes()
    for l in links:
        for h in range(len(hopt)):
            h0[h]=l
            if hopt[h]<l: h0[h]=hopt[h]
        N,F=sdcpf(N,h0=h0)
        N.save_nodes('Case_B_Link_'+str(l))
        save('./results/'+'Flows_Case_B_Link_'+str(l),F)

#scenario C
def Case_C(h0=None,betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]):
    if h0 == None: h0=get_quant(.99)
    N=Nodes()
    for b in betas:
        N,F=sdcpf(N,b=b,h0=h0)
        N.save_nodes('Case_C_Beta_'+str(b))
        save('./results/'+'Flows_Case_C_Beta_'+str(b),F)


#0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999
def Case_D(quants=[0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]):
    N=Nodes()
    for q in quants:
        h=get_quant(q)
#        if sum(h) >=1:
#            for i in range(len(h)):
#                if h[i]==0: h[i]=50
#            print h
        N,F=sdcpf(N,h0=h)
        N.save_nodes('Case_D_Quant_'+str(q))
        save('./results/'+'Flows_Case_D_Quant_'+str(q),F)

def biggestpair(H):
    H0=np.zeros((len(H))/2)
    for i in range(len(H0)):
        H0[i]=max(H[2*i],H[2*i+1])
    return H0

def Plot_A():
    betas=[0.0,0.25,0.5,0.75,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5]
    PlotA=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotA[j,0]=b
        N=Nodes(load_filename='Case_A_Beta_'+str(b)+'.npz')
        a=0
        d=0
        for i in N:
            a+=np.sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotA[j,1]=c
        j+=1
    save('./results/PlotA',PlotA)
    return PlotA


def Plot_B():
    links=np.arange(0.0,30000.1,1500.0)
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    hopt=get_quant(.99)
    h0=get_quant(.99)
    PlotB=np.zeros((len(links),2))
    j=0
    for l in links:
        for h in range(len(hopt)):
            h0[h]=l
            if hopt[h]<l: h0[h]=hopt[h]
        Hopt=biggestpair(h0)
        PlotB[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_B_Link_'+str(l)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotB[j,1]=c
        j+=1
    save('./results/PlotB',PlotB)
    return PlotB

def Plot_C():
    betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    Hop=get_quant(.99)
    Hopt=biggestpair(Hop)
    PlotC=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotC[j,0]=b*sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_C_Beta_'+str(b)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotC[j,1]=c
        j+=1
    save('./results/PlotC',PlotC)
    return PlotC

def Plot_D():
    quants=[0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    PlotD=np.zeros((len(quants),2))
    j=0
    for q in quants:
        Hop=get_quant(q)
        Hopt=biggestpair(Hop)
        PlotD[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_D_Quant_'+str(q)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotD[j,1]=c
        j+=1
    save('./results/PlotD',PlotD)
    return PlotD

def drawgrid(N):
    a,b,c,G=AtoKh_nonsparse(N)
    pos=nx.spring_layout(G)
    pos['FI']=[1.0,1.0]
    pos['SE']=[0.75,1.0]
    pos['NO']=[0.5,1.0]
    pos['DK']=[0.5,0.875]
    pos['PL']=[0.75,0.8]
    pos['GR']=[0.7,0.0]
    pos['BG']=[0.9,0.0]
    pos['SK']=[0.90,0.55]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RO']=[1.0,0.15]
    pos['RS']=[0.85,0.15]
    pos['BA']=[0.65,0.15]
    pos['HR']=[0.75,0.3]
    pos['IE']=[0.0,0.95]
    pos['GB']=[0.15,0.85]
    pos['FR']=[0.15,0.60]
    pos['ES']=[0.15,0.35]
    pos['PT']=[0.0,0.15]
    pos['BE']=[0.3,0.8]
    pos['NL']=[0.40,0.85]
    pos['LU']=[0.325,0.575]
    pos['DE']=[0.45,0.7]
    pos['CH']=[0.4,0.45]
    pos['IT']=[0.4,0.2]
    pos['AT']=[0.55,0.45]
    pos['SI']=[0.55,0.3]



    nx.draw_networkx_nodes(G,pos,node_size=1600,node_color='b',facecolor=(1,1,1))
    esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=500]
    emid=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>500 and d['weight']<=1500]
    elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500]

    nx.draw_networkx_edges(G,pos,edgelist=esmall,width=2)#1
    nx.draw_networkx_edges(G,pos,edgelist=emid,width=2)#3
    nx.draw_networkx_edges(G,pos,edgelist=elarge,width=2)#6
    nx.draw_networkx_labels(G,pos,font_size=24,font_color='w',font_family='sans-serif')
    axis('off')
    savefig("weighted_graph.png") # save as png
    show() # display

def plotbars(imports,exports): 
    order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]
    N=Nodes()
    names=['Ave.','','']
    for o in order:
        names.append(str(N[o].label))
    gcf().set_dpi(400)
    double_column_width = (2*3.425+0.236)*1.5 #extrawide!
    #gcf().set_size_inches([15,5])
    impcop=[0.38016,0,0]
    expcop=[0.38016,0,0]
    imp99=[0.369589,0,0]
    exp99=[0.369589,0,0]
    imp2x=[0.277789,0,0]
    exp2x=[0.277789,0,0]
    impact=[0.129063,0,0]
    expact=[0.129063,0,0]

    width=0.65
    for o in order:
        impcop.append(imports[o,0])
        imp99.append(imports[o,1])
        imp2x.append(imports[o,2])
        impact.append(imports[o,3])
        expcop.append(exports[o,0])
        exp99.append(exports[o,1])
        exp2x.append(exports[o,2])
        expact.append(exports[o,3])
    ind=np.arange(30)
    ax = subplot(111)
    rects1=bar(ind*1.25,expcop,width,align='center',color=(1,1,1))
    rects2=bar(ind*1.25+.1,exp99,width,align='center',color=(.8,.8,.8))
    rects3=bar(ind*1.25+.2,exp2x,width,align='center',color=(.6,.6,.6))
    rects4=bar(ind*1.25+.3,expact,width,align='center',color=(.4,.4,.4))
    pp = (rects1,rects2,rects3,rects4)
    pp_txtlabels = (r'Copper Plate',r'99th Quantile',r'Double of Actual',r'Actual Capacity')
    leg = legend(pp,pp_txtlabels,loc='upper right',ncol=len(pp));
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    gcf().set_size_inches([double_column_width*1.5,3*1.75])
    axis(ymin=0,ymax=0.75,xmin=amin(ind)-.875,xmax=amax(ind)+.875)
    xticks(ind*1.25+.65,names,rotation=60,ha='right',va='top')
    ylabel(r'Export Capacity Ratio')
    savefig('./figures/exports.pdf', dpi=300)
    #show()

def plot_allcases():
    plt.close()
    plota=load('./results/PlotA.npy')
#    plotb=load('./results/PlotB.npy')
#    plotc=load('./results/PlotC.npy')
    plotd=load('./results/PlotD.npy')
    plota[:,0]*=74.83
#    plotb[:,0]*=74.83
#    plotc[:,0]*=74.83
    plotd[:,0]*=74.83
    plota[:,1]*=100
    plotd[:,1]*=100
    ax=subplot(1,1,1)
    plt.axis([0, 900, 10, 30.1])
    plt.xlabel('Total Installed Transmission Capacity (GW)')
    plt.ylabel('Balancing Energy Required (% of Annual Consumption)')
    #plt.title(r'Reduction in required balancing energy with increased transmission capacity', fontsize=11)
    plt.text(-10,24.440991,'X',fontsize=20,color='r')
    plt.text(20,24.640991,'Individual Nodes')
    plt.text(839,14.9595,'X',fontsize=20,color='r')
    plt.text(74.83*9.20,16.31,'\'Copper Plate\'')
    plt.text(150,15.6,'Lower limit for Agreggated Nodes')
    p0=ax.plot([0,1500],[15.2595,15.2595],linestyle='dashed',color='k')
    xticks([])
    savefig('./figures/allcases1.pdf', dpi=300)
    plt.close()

    ax=subplot(1,1,1)
    plt.axis([0, 900.1, 10, 30.1])
    plt.xlabel('Total Installed Transmission Capacity (GW)')
    plt.ylabel('Balancing Energy Required (% of Annual Consumption)')
    #plt.title(r'Reduction in required balancing energy with increased transmission capacity', fontsize=11)
    plt.text(-10,24.440991,'X',fontsize=20,color='r')
    plt.text(20,24.640991,'Individual Nodes')
    plt.text(74.83*5.204590,15.2519634,'X',fontsize=20,color='r')
    plt.text(74.83*3.20-40,14.31,'99th Q -- Sufficiently Large Capacity')
    plt.text(839,014.9595,'X',fontsize=20,color='r')
    plt.text(74.83*9.20,16.31,'\'Copper Plate\'')
    p0=ax.plot([0,1500],[15.2595,15.2595],linestyle='dashed',color='k')
    savefig('./figures/allcases2.pdf', dpi=300)
    plt.close()

    ax=subplot(1,1,1)
    plt.axis([0, 900.1, 10, 30.1])
    plt.xlabel('Total Installed Transmission Capacity (GW)')
    plt.ylabel('Balancing Energy Required (% of Annual Consumption)')
    #plt.title(r'Reduction in required balancing energy with increased transmission capacity', fontsize=11)
    plt.text(-10,24.440991,'X',fontsize=20,color='r')
#    plt.text(20,.24640991,'Individual Nodes')
    plt.text(74.83*5.204590,15.2519634,'X',fontsize=20,color='r')
#    plt.text(74.83*3.20-40,.1631,'99th Q -- Sufficiently Large Capacity')
    plt.text(840,14.9595,'X',fontsize=20,color='r')
    p0=ax.plot([0,1500],[15.2595,15.2595],linestyle='dashed',color='k')
    ax.vlines(74.83,0.1,25,linestyle='dashed',color='k')
    ax.vlines(74.83*2,0.1,22,linestyle='dashed',color='k')
    ax.vlines(74.83*3,0.1,20,linestyle='dashed',color='k')
    ax.vlines(74.83*4,0.1,19,linestyle='dashed',color='k')
    ax.vlines(74.83*5,0.1,18.5,linestyle='dashed',color='k')
    plt.text(74.83*1,25.1,'x1 Actual Capacity')
    plt.text(74.83*2,22.1,'x2')
    plt.text(74.83*3,20.1,'x3')
    plt.text(74.83*4,19.1,'x4')
    plt.text(74.83*5,18.6,'x5')
    plote=np.append(plotd,[[900,15.2595]],0)
    p1,=ax.plot(plota[:,0],plota[:,1],label='Increase of Actual Caps.',linewidth=3)
    p4,=ax.plot(plote[:,0],plote[:,1],label='Reduction of 99% Q Caps.',linewidth=3)
    handles,labels=ax.get_legend_handles_labels()
    ax.legend(handles,labels)
    savefig('./figures/allcases3.pdf', dpi=300)
    plt.close()

    #show()

def find_balancing_reduction_quantiles(reduction=0.90,eps=1.e-3,guess=0.98):
    '''Loop over different quantile line capacities until the quantile
is found that leads to a reduction of balancing by <reduction>
percent, with a possible relative uncertainty of <eps>. Guess specifies
your first guess for the quantile.'''

    baltarget=.16197916
    olddist=0.
    balreal=0.
    N=Nodes()
    #N.set_alphas(0.7)
    #N.set_gammas(1.0)
    quant=guess # initial guess
    step=0.01
    print step
    while True:
        h=get_quant(quant)
        N,F=sdcpf(N,h0=h)
        a=0.; b=0.
        for i in N:
            a+=sum(i.balancing)
            b+=i.mean*i.nhours
        balreal=a/b
        reldist=abs(1.-balreal/baltarget)
        dist=baltarget-balreal
        if (reldist < eps):
            print '%15s %15s %15s %15s %15s' % ('distance','old distance','relative dist.','quantile','stepsize')
            print '%15.8f %15.8f %15.4f %15.4f %15.6f' % (dist, olddist, reldist,quant,step)
            del N
            break
        if (dist*olddist<0.): # sign change = we passed the perfect point! now reduce step size
            step=step/2.
        if dist<0:
            quant +=step
        if dist>0:
            quant -=step
        if (quant>=1.):
            step=step/2.
            quant=1.-step
        print '%15s %15s %15s %15s %15s' % ('distance','old distance','relative dist.','quantile','stepsize')
        print '%15.8f %15.8f %15.4f %15.4f %15.6f' % (dist, olddist, reldist,quant,step)
        olddist=dist
    del N
    return quant
# A after beta=15 is bad    
#Case_A()
#Case_B()
#Case_C(betas=[0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0])
#Case_D(quants=[0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89])
#Plot_A()
#Plot_B()
#Plot_C()
#Plot_D()
#plot_allcases()
