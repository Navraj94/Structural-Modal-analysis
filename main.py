# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 18:04:03 2020

@author: Naveen Raj
"""
import numpy as np
import vibs_analysis as va
import matplotlib.pyplot as plt
# =============================================================================

dof=3;
tinit=0;
step_size=.005;
tlimit=1000*step_size; 
tends=np.array([tinit,tlimit])
y0=np.zeros([dof,1])
yd0=np.zeros([dof,1])
y_int=np.vstack((y0,yd0))
alp=1/4
beta=.5
tspan=np.arange(tinit,tlimit,step_size)
#
# =============================================================================
#d,v,a=va.newmarkbeta_nonlin(dof,y_int,step_size,tends,alp,beta)
d1,v1,a1=va.newmarkbeta_lin(dof,y_int,step_size,tends,alp,beta)

omg,mod_shape,cita=va.SIMO_IBR_TDD(d1[:,100:600],tspan[100:600])