function rm=rhsMat(ff,r,N)
% in [0, 1]
rm=zeros(N+1,1);
[Pts,w]=JacobiGQ(0,0,20*N+1);
w=w/2;
Pts=Pts/2+1/2;
for j=0:N
    rm(j+1)=sum(ff(Pts.^r).*(Jac_Va(Pts,0,0,j,0)+Jac_Va(Pts,0,0,j+1,0)).*w);
end

end