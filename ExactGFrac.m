function gfrac=ExactGFrac(x,r,alpha,zn,wm)
%  Exact fractional derivative of order alpha(>0,<1) with respect to z(t)=t^zn and
%  weight w(t)=t^wn
% m = ceil(alpha);
% k=r+zn-1-m*zn;
% beta = alpha+1-m;
% if wn+r==0
%     gfrac=zeros(size(x));
%     return;
% else
%     C = gamma(wn+r+1)/gamma(wn+r+1-m)/zn^m*gamma((1+k)/zn)/gamma((1+k)/zn+1-beta);
% end
% gfrac = C*(x+eps).^(k+1-zn*beta);
gfrac =(x+eps) .^ (r - alpha * zn) * gamma((r + wm) / zn+1)  / gamma(-alpha+ 1 + (r + wm) / zn);