function J=Jac_Va(x,a,b,N,output_type)
if size(x,2)~=1
    x=x';
end
L=length(x);
J = ones(L,N+1);
if N==0
    J=J(:,1);
    return;
end
J(:,2)=(1+a)+(a+b+2)*(x/2-1/2);
if N==1
    if output_type==0
        J=J(:,2);
        return;
    else
        J=J(:,1:2);
        return;
    end
end
for n=2:N
    J(:,n+1) = ((2*n+a+b-1)*((2*n+a+b)*(2*n+a+b-2)*x+a^2-b^2).*J(:,n)-...
        2*(n+a-1)*(n+b-1)*(2*n+a+b)*J(:,n-1))/(2*n*(n+a+b)*(2*n+a+b-2));
end
if output_type==0
    J=J(:,N+1);
end
end