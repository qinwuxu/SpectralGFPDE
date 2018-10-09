function [rx,col_pts,M]=GFracMat(za,zb,N,alpha,w,invz)
% za=0;zb=1;
% N=22;
% alpha=2/3;
% w=inline('ones(size(x)).*x','x');
% invz=inline('x.^2','x');z=inline('x.^0.5','x');
r0=((JacobiGL(0,0,N)+1)/2);
col_pts=za+r0*(zb-za);  % -   - collocation points -  -  -
if alpha>0
FM=Jac_Frac_Diff_Shift(col_pts,0,0,N,alpha,za,zb,1);
else
    FM=(2/(zb-za))^(alpha)*Jac_Frac_Int(r0*2-1,0,0,N,-alpha,1);
    FM(1,:)=0;
end
%% choice of interpolation points 
r1 = r0;
interp_pts = r1*2-1;
Vand = Jac_Va(interp_pts,0,0,N,1);
%%
rx=invz(col_pts);
WL=diag(1./w(rx));
WL(1,1)=0;
WR=diag(w(rx));
M=WL*(FM/Vand)*WR;
% gfd=M*f(rx);
% plot(rx,gfd-ExactGFrac(rx,2,alpha,1/2,1),'.')
% function res=f(x)
% res=x.^2;
% end

end
