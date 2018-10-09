function jd=Jac_Der(x,a,b,N,m,output_type)
% x:   points, where derivatives should ne evaluated
% a,b: Indexes of jacobi polynomials
% N :  maximum order of Jacobi polnomials
% m :  order of  derivative ; 
% output_type: 'M'-output the matrix for all n=0,1,...N; 'V' - output a
%        vector for only n=N.

if size(x,2)~=1
    x=x';
end
J = Jac_Va(x,a+m,b+m,N-m,output_type);
if output_type == 'V' || output_type == 0
    jd = gamma(N+m+a+b+1)/2^m/gamma(N+a+b+1)*J;
else
    jd=zeros(length(x),N+1);
    for j=1+m:N+1
        jd(:,j)=gamma(j-1+m+a+b+1)/2^m/gamma(j-1+a+b+1)*J(:,j-m);
    end
end
    
