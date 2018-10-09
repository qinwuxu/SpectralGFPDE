function J = Jac_Frac_Int(x,a,b,N,alpha,output_type)
% Computing alpha-th order fractional integral of Jacobi polinomial of degree up to N
% a, b : J^{a,b}_n (x), n=0..N
% alpha : degree of derivative to be computed
% N : degree of Jacobi polinomial
% x : variables where derivatives need to be computed
% output_type: if output_type==0, only the result J_N will be outputed,
%   otherwise all the results J_n (n=0,..,N) will be outputed

if size(x,2)~=1
    x=x';
end
L=length(x);
J = zeros(L,N+1);
J(:,1) = (x+1).^alpha/gamma(alpha+1);
if N==0
    return;
end
J(:,2) = (a+b+2)/2*(x.*(x+1).^alpha/gamma(alpha+1)-alpha*(x+1).^(alpha+1)/gamma(alpha+2))+(a-b)/2*J(:,1);
if N==1
    return;
end
Jac_1 = Jac_Va(-1,a,b,N,1);
for n=2:N
    J(:,n+1) = (Aj(a,b,n-1)*x-Bj(a,b,n-1)-alpha*Aj(a,b,n-1)*Bhatj(a,b,n-1)).*J(:,n)-...
        (Cj(a,b,n-1)+alpha*Aj(a,b,n-1)*Ahatj(a,b,n-1)).*J(:,n-1)+...
        Aj(a,b,n-1)*(Ahatj(a,b,n-1)*Jac_1(n-1)+Bhatj(a,b,n-1)*Jac_1(n)+Chatj(a,b,n-1)*Jac_1(n+1))*...
        (x+1).^alpha/gamma(alpha);
    J(:,n+1) =J(:,n+1)/(1+alpha*Aj(a,b,n-1)*Chatj(a,b,n-1));
end
if output_type==0
    J=J(:,N+1);
end

    function aj = Aj(a,b,j)
        aj=(2*j+a+b+1).*(2*j+a+b+2)./(2*(j+1).*(j+a+b+1));
    end
    function bj = Bj(a,b,j)
        bj=(b.^2-a.^2).*(2*j+a+b+1)./(2*(j+1).*(j+a+b+1).*(2*j+a+b));
    end
    function cj = Cj(a,b,j)
        cj=(j+a).*(j+b).*(2*j+a+b+2)./((j+1).*(j+a+b+1).*(2*j+a+b));
    end
    function ahj = Ahatj(a,b,j)
        ahj = -2*(j+a).*(j+b)./((j+a+b).*(2*j+a+b).*(2*j+a+b+1));
    end
    function bhj = Bhatj(a,b,j)
        bhj = 2*(a-b)./((2*j+a+b).*(2*j+a+b+2));
    end
    function chj = Chatj(a,b,j)
        chj = 2*(j+a+b+1)./((2*j+a+b+1).*(2*j+a+b+2));
    end
end