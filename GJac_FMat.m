function M=GJac_FMat(beta,n,alp)
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


end