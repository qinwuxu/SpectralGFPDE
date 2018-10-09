function [x,u,ue]=PolySpectral(N,alp,lam)
% solving D^alp u(t) = lam* u(t) +1;
%         u(0) = 0;
GM=Jac_FMat(0,N,alp)*2^(alp-1);
DM=Jac_Inner_Mat(0,0,N)/2;
M=GM-lam*DM;
R=zeros(N+1,1);
R(1)=InnerJw(0,0,0)/2;
Uhat=M\R;

x=(0:.01:1)';
JM1=Jac_Va(2*x-1,0,0,N,1);
JM2=Jac_Va(2*x-1,0,0,N+1,1);
JM=JM1+JM2(:,2:end);
u=JM*Uhat;
ue=((x+0).^alp).*mlf(alp,alp+1,lam*(x+0).^alp,12);
