# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 18:27:06 2020
description: Functions needed for solution and analysis of 
             vibratory dynamical systems.
@author: Naveen Raj
"""
import numpy as np
import systmatrices as sm
    
def newmarkbeta_lin(dof,y_int,step_size,tends,alp,beta):

    tspan=np.arange(tends[0],tends[1],step_size)
    x=y_int[0:dof]
    dx=y_int[dof:2*dof+1]
    forc=np.zeros([dof,len(tspan)])
    disps=np.zeros([dof,len(tspan)])
    velos=np.zeros([dof,len(tspan)])
    accs=np.zeros([dof,len(tspan)])
    disps[:,0]=x.transpose()
    velos[:,0]=dx.transpose()
    m,k,damp,forc=sm.systmatrices(dof,tspan[0],y_int,step_size,tends) 
    accs[:,0]=np.matmul(np.linalg.inv(m),(forc[:,0]-np.matmul(damp,velos[:,0])-np.matmul(k,disps[:,0])))
    c1=1/(alp*step_size**2)
    c2=(1/(alp*step_size))
    c3=beta/(alp*step_size)
    c4=((1/(2*alp))-1)
    c5=beta/alp
    inver_mat=((c1*m)+(c3*damp)+k)
    matr=np.linalg.inv(inver_mat)
    for i in range(1,len(tspan[1:])):
        disps[:,i]=np.matmul(matr,(forc[:,i]+np.matmul(m,(c1*disps[:,i-1]+c2*velos[:,i-1]+c4*accs[:,i-1]))+np.matmul(damp,(c3*disps[:,i-1]+(c5-1)*velos[:,i-1]))+(c5-2)*(step_size/2)*accs[:,i-1]))
        accs[:,i]=c1*(disps[:,i]-disps[:,i-1])-c2*velos[:,i-1]-c4*accs[:,i-1]
        velos[:,i]=velos[:,i-1]+(1-beta)*step_size*accs[:,i-1]+beta*step_size*accs[:,i]
    return disps,velos,accs

# =============================================================================
#  WITH NEWTON-RAPHSON ITERATION FOR NON_LINEAR SYSTEMS
def newmarkbeta_nonlin(dof,y_int,step_size,tends,alp,beta,tol_rel=1e-6,maxit=20):

    tspan=np.arange(tends[0],tends[1],step_size)
    x=y_int[0:dof]
    dx=y_int[dof:2*dof+1]
    
    forc=np.zeros([dof,len(tspan)])
    disps=np.zeros([dof,len(tspan)])
    velos=np.zeros([dof,len(tspan)])
    accs=np.zeros([dof,len(tspan)])
    disps[:,0]=x.transpose()
    velos[:,0]=dx.transpose()
    m,k,damp,forc=sm.systmatrices(dof,tspan[0],y_int,step_size,tends) 
    accs[:,0]=np.matmul(np.linalg.inv(m),(forc[:,0]-np.matmul(damp,velos[:,0])-np.matmul(k,disps[:,0])))
    c1=1/(alp*step_size**2)
    c2=(1/(alp*step_size))
    #c3=beta/(alp*step_size)
    c4=((1/(2*alp))-1)
    #c5=beta/alp
    for i in range(1,len(tspan[1:])):
        m,k,damp,forc=sm.systmatrices(dof,tspan[i-1],np.vstack((disps[:,i-1],velos[:,i-1])),step_size,tends) 
        a1=(1/(alp*step_size**2)*m+(beta/(alp*step_size)*damp))
        a2=(1/(alp*step_size)*m+(beta/alp-1)*damp)
        a3=(1/(2*alp)-1)*m+(beta/(2*alp)-1)*damp
        #kcap=k+a1
        pcap=forc[:,i-1]+np.matmul(a1,disps[:,i-1])+np.matmul(a2,velos[:,i-1])+np.matmul(a3,accs[:,i-1])
        u=np.zeros([dof,maxit+1])
        j=0;
        u[:,j]=disps[:,i-1]
        fs=np.matmul(k,u[:,j])
        kt=np.copy(k)
        R=np.zeros([dof,maxit+2])
        R[:,j]=pcap-fs-np.matmul(a1,u[:,j])
        tol=2
        while tol > tol_rel and j<=maxit:
            kt=kt+a1
            delt=np.matmul(np.linalg.inv(kt),R[:,j])
            u[:,j+1]=u[:,j]+delt
            fs=np.matmul(k,u[:,j+1])
            R[:,j+1]=pcap-fs-np.matmul(a1,u[:,j+1])
            tol=np.linalg.norm(u[:,j+1]-u[:,j])
            j+=1
        disps[:,i]=u[:,j-1]
        velos[:,i]=(beta/(alp*step_size))*(disps[:,i]-disps[:,i-1])+(1-beta/alp)*velos[:,i-1]+step_size*(1-(beta/(2*alp)))*accs[:,i-1]
        accs[:,i]=((c1*(disps[:,i]-disps[:,i-1]))-(c2*velos[:,i-1])-(c4*accs[:,i-1]))
    return disps,velos,accs
# =============================================================================
# IBRAHIM TIME DOMAIN METHOD FOR MODAL ESTIMATION
def SIMO_IBR_TDD(disps,tspan):
    step_size=tspan[1]-tspan[0]
    X=disps[:,0:-3]
    Y=disps[:,1:-2]
    Z=disps[:,2:-1]
    V=np.vstack((X,Y))
    W=np.vstack((Y,Z))
    phi=np.dot(W,V.transpose())
    phi_hat=np.dot(V,V.transpose())
    phi_hat_i=np.linalg.inv(phi_hat)
    A=np.dot(phi,phi_hat_i)
    lamb,mod=np.linalg.eig(A)
    bet1=np.real(lamb)
    gam1=np.imag(lamb)
    a1=-1/(2*step_size)*np.log(bet1**2+gam1**2)
    b1=1/(step_size)*np.arctan(gam1/bet1)
    omg_ib1=np.sqrt(a1**2+b1**2);
    cita=a1/omg_ib1;
    return omg_ib1,mod,cita 

