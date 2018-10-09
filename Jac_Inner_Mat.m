function Lg=Jac_Inner_Mat(alp,beta,N)
Lg=zeros(N+1);

n=0;
IJ1=InnerJw(alp,beta,n);
IJ2=InnerJw(alp,beta,n+1);
Lg(n+1,n+1) = IJ1+IJ2;
Lg(n+2,n+1) = IJ2;

for n=1:N-1
    IJ1=InnerJw(alp,beta,n);
    IJ2=InnerJw(alp,beta,n+1);
    Lg(n+1,n+1) = IJ1+IJ2;
    Lg(n+0,n+1) = IJ1;
    Lg(n+2,n+1) = IJ2;
end

n=N;
IJ1=InnerJw(alp,beta,n);
IJ2=InnerJw(alp,beta,n+1);
Lg(n+1,n+1) = IJ1+IJ2;
Lg(n+0,n+1) = IJ1;
