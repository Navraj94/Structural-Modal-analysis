clc
clear all
close all
%%
% Input DOF file
% Define only the initial conditions and time steps in this file
% define the mass , damping , stiff in the systmatx.m
% things will take care by itself s_it back and relax
dof=3;
tinit=0;
step_size=.005;
tlimit=1000*step_size; % multiply by number of steps needed (12*step_size;) or directly input the end point
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
[m,c,ki,f1]=systmatx(t1(i),[disps2(:,i)';velo2(:,i)'],dof,tspan);    
accs2(:,i)=inv(m)*(f1-c*velo2(:,i)-ki*disps2(:,i));
end
D=eig(ki,m);
d_a=(sqrt(D));
%% =========================== PLOTTING FRF ===============================
k=1;
omg_r=1:0.01:120;
for i=1:length(omg_r)
    omg=omg_r(i);
    frf(:,:,k)=inv((-omg^2*m+1j*omg*c+ki));
    frf_11(k)=abs(frf(3,3,k));
    k=k+1;
end
% plot(omg_r,frf_11/norm(frf_11))
y_int=[y0 yd0];
tends=[tinit tlimit];
alpha=1/4;
beta=.5;
[disps1,velo1,accs1]=newmarkbeta(dof,y_int,alpha,beta,step_size,tends);

plot(disps1(1,:))
hold on
plot(disps2(1,:),'g-.')
err=disps2(:,:)-disps1(:,2:end);
%% ===================== IBRAHIM TIME DOMAIN METHOD =======================
X=disps1(:,101:597);
Y=disps1(:,102:598);
Z=disps1(:,103:599);
V=[X;Y];W=[Y;Z];
% [v1,d1]=eig(W*inv(V));
% lamb=cplxpair(diag(d1));  
ab=W*V';
ab2=V*V';
A=ab*inv(ab2);
[v2,d2]=eig(A);
lamb1=cplxpair(diag(d2));

% bet=real(lamb);
% gam=imag(lamb);
% a=1/(2*step_size)*log(bet.^2+gam.^2);
% b=1/(step_size)*atan(gam./bet);
% omg_ib=sqrt(a.^2+b.^2);
% cita=a./omg_ib;


bet1=real(lamb1);
gam1=imag(lamb1);
a1=-1/(2*step_size)*log(bet1.^2+gam1.^2);
b1=1/(step_size)*atan(gam1./bet1);
omg_ib1=sqrt(a1.^2+b1.^2)
cita=a1./omg_ib1;

% [omg_ibR1 omg_ib]
%% ========================= NEWMARK_BETA =================================
% y_int=[y0 yd0];
% tends=[tinit tlimit];
% alp=1/4;
% beta=.5;
% [disps1,velo1,accs1]=newmarkbeta(dof,y_int,alp,beta,step_size,tends);
% 

% y_int=[y0 yd0];
% tends=[tinit tlimit];
% alp=1/4;
% beta=.5;
% [disps3,velo3,accs3]=newmarkbeta_2(dof,y_int,alp,beta,step_size,tends);