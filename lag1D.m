function lag=lag1D(x,N,zx)
if nargin<3
    zx = JacobiGQ(0,0,N)';
end
[m,n]=size(x);
if m==1&&n>1
    x=x';
end
[m,n]=size(zx);
if m>1&&n==1
    zx=zx';
end
Lx=length(x);
lag=zeros(Lx,N+1);
for i=1:N+1
    deno = prod(zx(i)-[zx(1:i-1),zx(i+1:N+1)]);
    lag(:,i)=prod(x*ones(1,N)-ones(Lx,1)*[zx(1:i-1),zx(i+1:N+1)],2)/deno;
end