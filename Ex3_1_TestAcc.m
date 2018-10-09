za=0;zb=1;
alpha=0.8;
w=inline('ones(size(x))','x');
invz=inline('x.^(2/3)','x');z=inline('x.^(3/2)','x');
zn=3/2;
wn=0;
n=[4:12];
p=15;q=18;
EA=[];

for alpha=[0.2 0.5 0.8]
    e1=[];e2=[];
for j=1:length(n)
    N=n(j);
    [rx1,r1,M1]=GFracMat(za,zb,N,alpha,w,invz);
%     [rx2,r2,M2]=TestGFracMat(za,zb,N,alpha,w,z,invz);
    y1= rx1.^3/prod(1:2)       -  rx1.^6/prod(1:4) ...
       +rx1.^9/prod(1:6)       -  rx1.^12/prod(1:8) ...
       +rx1.^15/prod(1:10) ;%    -   rx1.^18/prod(1:12);%-rx1.^q;
%     y2=rx2.^p-rx2.^q;
    % edy1=gamma(7)/gamma(7-alpha)*rx1.^(3-alpha/2);%
    edy1=  ExactGFrac(rx1,3,alpha,zn,wn)/prod(1:2)   -  ExactGFrac(rx1,6,alpha,zn,wn)/prod(1:4) ...
         + ExactGFrac(rx1,9,alpha,zn,wn)/prod(1:6) -  ExactGFrac(rx1,12,alpha,zn,wn)/prod(1:8) ...
         + ExactGFrac(rx1,15,alpha,zn,wn)/prod(1:10);% -  ExactGFrac(rx1,18,alpha,zn,wn)/prod(1:12);
%     edy2=ExactGFrac(rx2,p,alpha,zn,wn)-ExactGFrac(rx2,q,alpha,zn,wn);
    e1=[e1,max(2^4*abs(edy1-M1*y1))];
%     e2=[e2,max(abs(edy2-M2*y2))];
end
EA=[EA;e1];
end
% plot(rx1,edy1,'.',rx1,M1*y1,'r')
semilogy(n,EA(1,:),'s-',n,EA(2,:),'^-',n,EA(3,:),'*-','LineWidth',2)
legend('alpha=0.2','alpha=0.5','alpha=0.8')
xlabel('N');
ylabel('maximum error')