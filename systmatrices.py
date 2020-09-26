# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 19:11:25 2020

@author: Naveen Raj
"""

def systmatrices(dof,t,y1,step_size,tends):
    import numpy as np
    tspan=np.arange(tends[0],tends[1],step_size)
    x=y1[0:dof,:]
    dx=y1[dof+1:2*dof+1,:]
    # =============================================
    m=np.zeros([dof,dof])
    m=np.array([[100,0,0],[0,10,0],[0,0,10]])
    # =============================================
    c=np.zeros([dof,dof])
    c=np.array([[4,-2,0],[-2,4,-2],[0,-2,2]])*100
    # =============================================
    k=np.zeros([dof,dof])
    k=np.array([[8,-4,0],[-4,8,-4],[0,-4,4]])*1e4
    # =============================================
    f=np.zeros([dof,len(tspan)])
    for i in range(1000):
        if tspan[i] <=0.5:
            f[1,i]=1000
    
    Mass,K_stiff,C_damp,Forcs=m,k,c,f
    return Mass,K_stiff,C_damp,Forcs