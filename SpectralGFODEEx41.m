% Example:
% Solving the following ordinary fractional differential equation
%   D_{0,t}^\alpha u(t) = lambda*u(x) + f(t)
function [rx,M,U,exact]=SpectralGFODEEx41(N,alpha,zn,wn,p)
%zn=1/3;
%wn=.7;p=6;
% N=20;alpha=1;
[rx,r,D]=GFracMat(0,1,N,alpha,@w,@invz);
M=D-diag(lambda(rx));
M(1,:)=0;M(1,1)=1;
rhs=F(rx,alpha);
rhs(1)=0;
U=M\rhs;

exact = rx.^p;%AprSin(rx,20);%gamma(3)/gamma(3+alpha)*(r+1).^(2+alpha);

% plot(rx,U-exact,'.')
    function f=F(x,alpha)
        %     f=GFAprSin(x,alpha,zn,wn,20);%
            f=ExactGFrac(x,p,alpha,zn,wn)-lambda(x).*x.^p;
%     f=0*sin(2*pi*x.^alpha);
%         f=gamma(8.05)/gamma(8.05-alpha)*x.^(3.7-alpha*2/3).*(1+x);
    end
    function res=w(x)
        res=x.^wn;
    end
    function res=invz(x)
        res=x.^(1/zn);
    end
    function lam=lambda(x)
        lam=(1+x);%gamma(8.05)/gamma(8.05-alpha)*x.^(1-alpha*2/3);
    end
end
