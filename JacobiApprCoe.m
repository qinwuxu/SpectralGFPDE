function C=JacobiApprCoe(f,a,b,N)

[x,w]=lglnodes(2*N);
C=zeros(N+1,1);
J=Jac_Va(x,a,b,N,1);
for j=1:N+1
    deno=sum(J(:,j).^2.*w);
    C(j) = sum(f(x).*J(:,j).*w)/deno;
end
