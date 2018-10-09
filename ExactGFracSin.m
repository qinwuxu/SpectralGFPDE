function gsin=ExactGFracSin(rx,n,alpha,zn,wn)

gsin=0;

for j=1:2:n
    gsin=gsin-(-1)^(j/2+1/2)*ExactGFrac(rx,j,alpha,zn,wn)/prod(1:n);
end