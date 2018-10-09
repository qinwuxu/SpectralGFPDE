function [X,t,U]=SpaceTimeSpectral(M,N)

% L=(M+1)*(N+1);
x=JacobiGL(0,0,M);
t=JacobiGL(0,0,N)/2+1/2;
[X,T]=meshgrid(x,t);
Mx=Jac_Frac_Diff(x,0,0,M,2,1);
Vt=Jac_Va(2*t-1,0,0,N,1);
Mj=Jac_Frac_Int(2*t-1,0,0,N,1,1)/2;
Mt=Mj/Vt;

% MM2=kron(Mx,eye(N+1));
% MMb=kron(eye(M+1),Mt);
% MA=MMb*MM2;
% MA=(MMb-MM2);
%% Boundary and initial condition
% MA(1:N+1,:)=0; % for left boundary condition
% MA(1:N+1,1:N+1)=eye(1+N);
% MA(1:N+1:L,:)=0;  %for initial condition 
% MA(1:N+1:L,1:N+1:L)=eye(1+M);
% MA(L-N:L,:)=0; % for right boundary
% MA(L-N:L,L-N:L)=eye(N+1);
% U0=zeros(L,1);   %
% U0(1:N+1:L)=u0(rx); % initial condition
 %%
% U0=zeros(L,1);
% U0(1:M+1)=initU0(x);
% U=(eye(L)-MA)\U0;
% U=reshape(U,N+1,M+1);
U0=zeros(N+1,1);
U0(:)=1;
U=(eye(N+1)-Mt)\U0;

    function u0=initU0(x)
        u0=1-x.^2;
    end
end