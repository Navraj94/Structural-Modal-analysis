function [disps,velo,accs]=newmarkbeta(dof,y0,alp,beta,step_size,tends)
%           1.first coulmn of all the matrices correspond to the initial
%             values(in case of non zero initail values pls write another
%             step to assign them to the first column)
%
%           2.refer algorithm in ss.rao(mechanical vibrations)


% [row col]=size(dof);
forc=zeros(dof,1);
tspan=tends(1):step_size:tends(2);
% int_disp=zeros(length(k),1);
% int_veloct=zeros(length(k),1);
disps_t=zeros(dof,length(tspan));
accs_t=zeros(dof,length(tspan));
velo_t=zeros(dof,length(tspan));
disps_t(:,1)=y0(:,1);
velo_t(:,1)=y0(:,2);
% velo_t(1:3:76,1)=V_o;
% disps_t(3:3:78,1)=intcond(:,1);
%  disps_t(2:3:77,1)=-tand(10);
% disps_t(77,1)=tand(kkk);
[m,damp,k,f1]=systmatx(tends(1),[disps_t(:,1);velo_t(:,1)],dof,tspan);
accs_t(:,1)=inv(m)*(f1-(damp*velo_t(:,1))-(k*disps_t(:,1)));
for i=1:length(tspan)
    [m,damp,k,forc(:,i)]=systmatx(tspan(i),[disps_t(:,i);velo_t(:,i)],dof,tspan);
c1=1/(alp*step_size^2);
c2=(1/(alp*step_size));
c3=beta/(alp*step_size);
c4=((1/(2*alp))-1);
c5=beta/alp;

inver_mat=((c1*m)+(c3*damp)+k);
matr=inv(inver_mat);

    
%     hold on
%     plot(tspan(i),forc(77,i),'r*')
    
    %   if i == 1
    disps_t(:,i+1)=(matr*(forc(:,i)+(m*((c1*disps_t(:,i))+(c2*velo_t(:,i)+c4*(accs_t(:,i)))))+(damp*((c3*disps_t(:,i))+((c5-1)*velo_t(:,i))+(c5-2)*(step_size/2)*(accs_t(:,i))))));
    %     disps_t(:,i+1)=(matr*{forc(:,i)+(m*((c1*int_disp)+(c2*int_velo)+c4*(int_accs)))+(damp*((c3*int_disp)+((c5-1)*int_velo))+(c5-2)*(steps/2)*(int_accs))});
    accs_t(:,i+1)=((c1*(disps_t(:,i+1)-disps_t(:,i)))-(c2*velo_t(:,i))-(c4*accs_t(:,i)));
    velo_t(:,i+1)=(velo_t(:,i)+((1-beta)*step_size*accs_t(:,i))+(beta*step_size*accs_t(:,i+1)));
    
end
% figure
%     hold on
%     plot(fact)
disps=disps_t;
velo=velo_t;
accs=accs_t;

 %% rectangular
    
    %
%     %% Righ triangular pulse
%     switch count
%         case 1
%             if tspan(i) <= tlimit2
%                 forc(:,i)=tspan(i)/tlimit2*forceda;
%             else
%                 forc(:,i)=zeros(row,1);
%             end
%             %
%             % triangular pulse
%         case 2
%             if tspan(i) <= tlimit2
%                 
%                 if tspan(i) <= tlimit/2
%                     forc(:,i)=tspan(i)/tlimit2*2*forceda;
%                 else
%                     forc(:,i)=(1-tspan(i)/tlimit2)*2*forceda;
%                 end
%             else
%                 forc(:,i)=zeros(row,1);
%             end
%         otherwise
%             if tspan(i) <= tlimit2
%                 forc(:,i)=forceda;
%             else
%                 forc(:,i)=zeros(row,1);
%             end
%     end