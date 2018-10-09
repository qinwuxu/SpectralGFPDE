function [X,T,U]=Stretch_Diffusion_Int(t,alp,beta,M,N)
w2=inline('ones(size(x))','x');
z2=inline('x','x');
L=(M+1)*(N+1);
x=JacobiGL(0,0,M);

[rx,~,M2]=GFracMat(-1,1,M,2,w2,z2);
[rt,~,Mb]=GFracMat(z(0),z(t),N,-beta,w2,@invz);
[X,T]=meshgrid(rx,rt);

MM2=kron(M2,eye(N+1));
MMb=kron(eye(M+1),Mb);
MA=MMb*MM2;

%% Boundary condition --  no initial is needed
MA(1:N+1,:)=0; % for left boundary condition
MA(L-N:L,:)=0; % for right boundary

U0=reshape(repmat(u0(x'),N+1,1),[],1);
U0(1:N+1)=0;
U0(L-N:L,:)=0;

U=(eye(L)-MA)\U0;
U=reshape(U,N+1,M+1);

    function res=u0(x)
%         res=1-x.^4;%
        res=exp(-10*x.^2)-exp(-10);
    end

    function res=z(x)
        res=x.^(beta/alp);
    end
    function iz=invz(x)
        iz=x.^(alp/beta);
    end

end

