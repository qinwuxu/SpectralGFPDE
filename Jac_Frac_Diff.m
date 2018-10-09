function JFD = Jac_Frac_Diff(x,a,b,N,alpha,output_type)
% Computing alpha-th order fractional derivative of Jacobi polinomial of degree n
% a, b : J^{a,b}_n (x)
% alpha : degree of derivative to be computed
% n : degree of Jacobi polinomial
% x : variables where derivatives need to be computed
% output_type: if output_type==0, only the result J_N will be outputed,
%   otherwise all the results J_n (n=0,..,N) will be outputed

if size(x,2)~=1
    x=x';
end
L=length(x);

k=ceil(alpha);
d = gamma((0:N)+k+a+b+1)./(2^k*gamma((0:N)+a+b+1));

if output_type==0
    if N<k
        JFD = zeros(size(x));
        return;
    end
    JFD = d(N+1)*Jac_Frac_Int(x,a+k,b+k,N-k,k-alpha,0);
else
    JFD = zeros(L,N+1);
    JFD(:,k+1:N+1) = Jac_Frac_Int(x,a+k,b+k,N-k,k-alpha,1);
    for j=k+1:N+1
        JFD(:,j)=d(j)*JFD(:,j);
    end
end
end