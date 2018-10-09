
alp=1;beta=1;
[X,T,U]=Stretch_Diffusion_Int(0.5,alp,beta,50,50);
figure;surf(X,T,U,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp
set(gca,'FontSize',12);
xlabel('x');ylabel('t');zlabel('u(x,t)')

alp=1.5;beta=1;
[X,T,U]=Stretch_Diffusion_Int(0.5,alp,beta,50,50);
figure;surf(X,T,U,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp
set(gca,'FontSize',12);
xlabel('x');ylabel('t');zlabel('u(x,t)')

alp=0.5;beta=1;
[X,T,U]=Stretch_Diffusion_Int(0.5,alp,beta,50,50);
figure;surf(X,T,U,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp
set(gca,'FontSize',12);
xlabel('x');ylabel('t');zlabel('u(x,t)')

alp=0.6;beta=0.6;
[X,T,U]=Stretch_Diffusion_Int(0.5,alp,beta,50,50);
figure;surf(X,T,U,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
shading interp
set(gca,'FontSize',12);
xlabel('x');ylabel('t');zlabel('u(x,t)')