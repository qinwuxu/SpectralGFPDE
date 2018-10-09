function [r,U,exact] = SpectralFracODE_Shift(N,alpha)

Mx=pi;
r = Mx*(JacobiGL(0,0,N)+1)/2;
FJ = Jac_Frac_Diff_Shift(r,0,0,N,alpha,0,Mx,1);
Vand = Jac_Va_Shift(r,0,0,N,0,Mx,1);
M= FJ-Vand;
M(1,:)=Vand(1,:);
rhs = F(r,alpha);
rhs(1) = 0;
U = M\rhs;
U=Vand*U;
exact = Sin_Tay(r,30);%gamma(3)/gamma(3+alpha)*(r+1).^(2+alpha);

% plot(r,U-exact)
function f=F(t,alpha)
    f=Sin_Frac_Diff(t,alpha,30)-Sin_Tay(t,30);
end
end