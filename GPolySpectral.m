function [x,u,ue]=GPolySpectral(N,alp,lam)
% solving D^alp u(t) = lam* u(t) +1;
%         u(0) = 0;
FM=GJac_FracMat(N,alp)/2^(1-alp);
JM=GJac_StiffMat(N,alp)/2;
M=FM-lam*JM;
R=zeros(N+1,1);
R(1)=InnerJw(0,0,0)/2;
Uhat=M\R;

x=(0:.01:1)';
JM1=repmat((2*x).^alp,1,N+1).*Jac_Va(2*x-1,-alp,alp,N,1);
u=JM1*Uhat;
ue=((x+0).^alp).*mlf(alp,alp+1,lam*(x+0).^alp,12);


    function G=GJac_FracMat(n,alp)
        G=zeros(n+1);
        for j=0:n
            jj=1:j;
            G(j+1,j+1)=prod((jj+alp)./jj)*gamma(1+alp)*InnerJw(0,0,j);
        end
    end
    function M=GJac_StiffMat(n,alp)
        M=zeros(n+1);
        [Pts,w]=JacobiGQ(0,0,20*n+1);
        for m=0:n
            for j=0:n
                M(m+1,j+1)=sum(Jac_Va(Pts,0,0,m,0).*...
                    (1+Pts).^alp.*Jac_Va(Pts,-alp,alp,j,0).*w);
            end
        end
    end

end