function [x,u,ue]=EK_Spectral(N,alp,r,lam)
% D^\alpha u = lam* u + 1
GM=Jac_GMat(0,0,N,alp,r)/2^(1-alp*r);
DM=Jac_Inner_Mat(0,0,N)/2;
M=GM-lam*DM;
R=zeros(N+1,1);
R(1)=InnerJw(0,0,0)/2;
Uhat=M\R;

x=(0:.01:1)';
zx=x.^(1/r);
JM1=Jac_Va(2*zx-1,0,0,N,1);
JM2=Jac_Va(2*zx-1,0,0,N+1,1);
JM=JM1+JM2(:,2:end);
u=JM*Uhat;
ue=((x+0).^alp).*mlf(alp,alp+1,lam*(x+0).^alp,12);
% e=u-ue';
% plot(x,u,x,ue','.')

