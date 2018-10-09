function M=Jac_FMat(beta,n,alp)
% compute inner product :
%   M_{k,n}=( D^alp J^{beta,0}_n(x), J^{beta,0}_k(x) )
M=zeros(n+1);
[Pts,w]=JacobiGQ(beta,0,20*n+1);
for m=0:n
    for j=0:n
        M(m+1,j+1)=sum((Jac_Va(Pts,beta,0,m,0)+Jac_Va(Pts,beta,0,m+1,0)).*...
            (Jac_Frac_Diff(Pts,beta,0,j,alp,0)+Jac_Frac_Diff(Pts,beta,0,j+1,alp,0)).*w);
    end
end
% for m=0:k
%     for j=ceil(alp):n
%         ele=ele+(-1)^(k+m+j+n)*2^(1-alp)*GOG(m,k+beta+1)/(m+j-alp+1)*...
%             GOG(k-m,m+1)*GOG(j,n+beta+1)*GOG(n-j,j-alp+1)*Gam2(n,alp);
%     end
% end
% 
%     function g=GOG(m,a)
%         % computing gamma(m+a)/gamma(a)/gamma(m+1)
% %         if m==0
% %             g=1;
% %         else
%             mm=1:m;
%             g=prod((a+mm-1)./mm);
% %         end
%         
%     end
%     function g2=Gam2(a,alp)
%         % computing Gamma(a+1)/Gamma(a+1-alp)
%         inta=floor(a);
%         pa=a-inta;
%         sinta=1:inta;
%         g2=prod((sinta+pa)./(sinta+pa-alp))*gamma(1+pa)/gamma(1+pa-alp);
%     end
end