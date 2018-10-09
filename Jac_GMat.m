function M=Jac_GMat(beta,gamm,N,alp,r)
%

M=zeros(N+1);

for n=0:N
    for k=0:N
        M(k+1,n+1)=Jac_GMat_Ele(beta,gamm,n,alp,r,k)+...
        Jac_GMat_Ele(beta,gamm,n,alp,r,k+1)+...
        Jac_GMat_Ele(beta,gamm,n+1,alp,r,k)+...
        Jac_GMat_Ele(beta,gamm,n+1,alp,r,k+1);
    end
end