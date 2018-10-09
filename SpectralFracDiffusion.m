function [r,U,exact]=SpectralFracDiffusion(N,alpha,T)
r = JacobiGL(0,0,N)/2+1/2;
CFL=0.0001;
dt = CFL*(r(2)-r(1))^alpha;
Nstep=floor(T/dt);
dt = T/Nstep;
FJ = Jac_Frac_Diff_Shift(r,0,0,N,alpha,0,1,1);
Vand = Jac_Va_Shift(r,0,0,N,0,1,1);
M = FJ/Vand;
U=exactU(r,0,alpha);
resu=0;
rk4a=0;rk4b=0;rk4c=0;
RK4_Coe;x=r;
fxt=x.^(3 )+gamma(4 )/gamma(4-alpha)*x.^(3-alpha);
for nt=0:Nstep
    t=nt*dt;
    for it=1:5
        tlocal=t+rk4c(it)*dt;
        rhside=M*U+F(r,tlocal,alpha,fxt);
        resu=rk4a(it)*resu+dt*rhside;
        U=U+rk4b(it)*resu;
        U([1 N+1])=exactU(r([1 N+1]),t+rk4c(it+1)*dt,alpha);
    end
end
exact = exactU(r,T,alpha);


    function f=F(x,t,alpha,fxt)
        f=-exp(-t)*(fxt);%-exp(-t)*Sin_Tay(x+1,30)-exp(-t)*Sin_Frac_Diff(x+1,alpha,30);
    end
    function U=exactU(x,t,alpha)
        U=exp(-t)*x.^(3);%exp(-t)*Sin_Tay(x+1,30);
    end
end