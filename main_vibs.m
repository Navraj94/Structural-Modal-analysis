clc
clear all
close all
%%
% Input DOF file
% Define only the initial conditions and time steps in this file
% define the mass , damping , stiff in the systmatx.m
% things will take care by itself sit back and relax
dof=3;
tinit=0;
step_size=.01;
tlimit=500*step_size; % multiply by number of steps needed (12*step_size;) or directly input the end point
tspan=[0:step_size:tlimit];
y0=zeros(dof,1);
yd0=zeros(dof,1);
y_int=[y0;yd0];
tends=[tinit tlimit];
[t1,y1]=ode45(@(t,y)odefuns_vibr(t,y,dof,tspan),tspan,y_int);

disps2=y1(:,1:dof);
velo2=y1(:,dof+1:end);
disps2=disps2';
velo2=velo2';
for i=1:length(t1)
[m,c,k,f1]=systmatx(t1(i),[disps2(:,i)';velo2(:,i)'],dof,tspan);    
accs2(:,i)=inv(m)*(f1-c*velo2(:,i)-k*disps2(:,i));
end


%% ================== NEWMARK_BETA ====================
y_int=[y0 yd0];
tends=[tinit tlimit];
alpha=1/4;
beta=.5;
[disps1,velo1,accs1]=newmarkbeta(dof,y_int,alpha,beta,step_size,tends);


y_int=[y0 yd0];
tends=[tinit tlimit];
alpha=1/4;
beta=.5;
[disps3,velo3,accs3]=newmarkbeta_2(dof,y_int,alpha,beta,step_size,tends);