function [mass,damps,stiffs,forces]=systmatx(t,y1,dof,tspan)
x=y1(1:dof);
dx=y1(dof+1:end);
x=x';  % current time step disp values
dx=dx'; % current time step velo values
%%  ===================== Mass matrix formulations =======================
m=zeros(dof);
m=[100 0 0; 0 10 0;0 0 10];

%% ====================== Stiffness matrix formulations ==================
k=zeros(dof);
k=[8 -4 0; -4 8 -4;0 -4 4]*1e4;
%% ====================== Damping matrix formulations ====================
% c=zeros(dof);
% [v,d]=eig(m,k);
% omg=sqrt(diag(d));
% c1=0.01;c2=0.005;
% A=[1/omg(1) omg(1);1/omg(2) omg(2)];
% B=[2*c1;2*c2];
% alpB=inv(A)*B;
% c=alpB(1)*m+alpB(2)*k;

c=100*[4 -2 0;-2 4 -2;0 -2 2];

%% ====================== Force vector formulations ======================
% forc=zeros(dof,length(tspan));
% forc(2,:)=10;
if t< 0.5
f1=zeros(dof,1);
f1(2)=1000;
else
    f1=zeros(dof,1);
end
% for i=1:dof(1) % interpolated force vector value
%     f1(i,1)=forc(i,1)*sin(t); 
% end
% for i=1:length(tspan) % turn on to find accelerations 
% accs2(:,i)=inv(m)*(f1-c*dx(:,i)-k*x(:,i));
% end
mass=m;
damps=c;
stiffs=k;
forces=f1;
