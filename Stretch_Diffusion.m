function [X,T,U]=Stretch_Diffusion(t,alp,beta,M,N)
w2=inline('ones(size(x))','x');
z2=inline('x','x');
L=(M+1)*(N+1);
x=JacobiGL(0,0,M);

[rx,~,M2]=GFracMat(-1,1,M,2,w2,z2);
[rt,~,Mb]=GFracMat(z(0),z(t),N,beta,w2,@invz);
[X,T]=meshgrid(rx,rt);

MM2=kron(M2,eye(N+1));
MMb=kron(eye(M+1),Mb);
% MA=MMb*MM2;
MA=(MMb-MM2);
%% Boundary and initial condition
MA(1:N+1,:)=0; % for left boundary condition
MA(1:N+1,1:N+1)=eye(1+N);
MA(1:N+1:L,:)=0;  %for initial condition 
MA(1:N+1:L,1:N+1:L)=eye(1+M);
MA(L-N:L,:)=0; % for right boundary
MA(L-N:L,L-N:L)=eye(N+1);
U0=zeros(L,1);   %
U0(1:N+1:L)=u0(rx); % initial condition
 %%
U=(MA)\U0;
U=reshape(U,N+1,M+1);

    function res=u0(x)
%         res=exp(-100*x.^2);
        res=1-x.^2;
    end

    function res=z(x)
        res=x.^(beta/alp);
    end
    function iz=invz(x)
        iz=x.^(alp/beta);
    end

end

