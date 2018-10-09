function fr=J0_Frac(x,a,b,N,alp,RLorC,LorR,output_type)
% fr=J0_Frac(x,a,b,N,alp,RLorC,LorR,output_type)
% computing fractional derivative(alp>0) or integral(alp<0) of Jacobi polynomials
% x:   points, where fractional derivative should ne evaluated
% a,b: Indexes of jacobi polynomials
% N :  maximum order of Jacobi polnomials
% alp : order of fractional derivative (alp >0); If alp < 0, compute
%        fractional integral of order |alp|.
% RLorC : Type of fractional derivative. 'C'-Caputo; 'R':Riemann-Liouville
% output_type: 'M'-output the matrix for all n=0,1,...N; 'V' - output a
%        vector for only n=N.

if nargin < 5
    error('Not enough input arguments!');
elseif nargin < 6
    RLorC = 'C';
    LorR = 'L';
    output_type = 'M';
elseif nargin < 7
    LorR = 'L';
    output_type = 'M';
elseif nargin < 8
    output_type = 'M';
end
if size(x,2)~=1
    x=x';
end
if LorR == 'R'
    X=(1-x).^(-alp);
    J=Jac_Va(x,-alp,b+alp,N,1);
else
    X=(1+x).^(-alp);
    J=Jac_Va(x,a+alp,-alp,N,1);
end

fr = (X*(factorial(0:N)./gamma((0:N)+1-alp))).*J;

if sum(RLorC=='C') && alp >0
    falp=floor(alp);
    if LorR == 'R'
        for j=0:falp
            JD=Jac_Der(1,0,b,N,j,'M');
            fr = fr-(-1)^j/gamma(1+j-alp)*(1-x).^(j-alp)*JD;
        end
    else
        for j=0:falp
            JD=Jac_Der(-1,a,0,N,j,'M');
            fr = fr-1/gamma(1+j-alp)*(1+x).^(j-alp)*JD;
        end
    end
end
if output_type == 'V'
    fr=fr(:,N+1);
else
    fr(1,:) = 0;
end
    
