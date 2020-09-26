function dydt=odefuns_vibr(t,y,dof,tspan)
[m,c,k,f1]=systmatx(t,y,dof,tspan);
A=[zeros(dof) eye(dof);-inv(m)*k -inv(m)*c];
F=[zeros(dof,1);inv(m)*f1];
dydt=A*y+F;