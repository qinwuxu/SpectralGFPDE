function dsin=GFAprSin(x,alpha,zn,wn,N)
dsin=zeros(size(x));
for k=1:2:N
%     if k<alpha
%         xk=0;
%     else
%         xk=x.^(k-alpha);
%     end
    dsin=dsin+(-1)^((k-1)/2)/prod(1:k)*ExactGFrac(x,k,alpha,zn,wn);
end