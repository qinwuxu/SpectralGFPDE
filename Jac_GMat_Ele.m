function ele=Jac_GMat_Ele(beta,gamm,n,alp,r,k)
% compute inner product :
%   ( D^alp_z (J^{beta,gamm}_n(x)+J_{n+1}), (J^{beta,gamm}_k(x)+J_{k+1}) )

ele=0;        

for m=0:k
    for j=ceil(alp):n
        ele=ele+(-1)^(k+m+j+n)*2^(1-alp*r)*GOG(m,k+beta+gamm+1)/(m+j-alp*r+1)*...
            GOG(k-m,m+gamm+1)*GOG(j,n+beta+gamm+1)*GOG(n-j,j+gamm+1)*Gam2(j/r,alp);
    end
end

    function g=GOG(m,a)
        % computing gamma(m+a)/gamma(a)/gamma(m+1)
%         if m==0
%             g=1;
%         else
            mm=1:m;
            g=prod((a+mm-1)./mm);
%         end
        
    end
    function g2=Gam2(a,alp)
        % computing Gamma(a+1)/Gamma(a+1-alp)
        inta=floor(a);
        pa=a-inta;
        sinta=1:inta;
        g2=prod((sinta+pa)./(sinta+pa-alp))*gamma(1+pa)/gamma(1+pa-alp);
    end
end