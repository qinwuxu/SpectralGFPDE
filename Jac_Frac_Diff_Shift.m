function FJ = Jac_Frac_Diff_Shift(x,a,b,N,alpha,min,Max,output_type)
L=Max-min;
x=2*(x-min)/L-1;
FJ = (2/L)^alpha*Jac_Frac_Diff(x,a,b,N,alpha,output_type);
