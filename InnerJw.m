function res=InnerJw(alp,beta,n)
nn=1:n;
res=2^(alp+beta+1)/(2*n+alp+beta+1)*prod((nn+alp)./nn)*gamma(alp+1)*...
    prod((nn+beta)./(nn+alp+beta))*gamma(beta+1)/gamma(alp+beta+1);
end