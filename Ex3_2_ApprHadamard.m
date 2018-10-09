function gfd=ApprHadamard;%(ff,invz,ww,alp,x,N)
mu=0;
ff=inline('sin(pi*x)','x');
invz=inline('exp(x)','x');
x=1:.02:2;
ww=inline(strcat('x.^',num2str(mu)),'x');

%% ---- 1st
mu=0;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-1;
%% reference solution
N=50;
Coe=JacobiApprCoe(@gg,0,0,N);
Lx=(max(log(x))-min(log(x)));
X=2*(log(x)-min(log(x)))/Lx-1;
J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
gfd0=zeros(size(x'));
for j=0:N
    gfd0=gfd0+Coe(j+1)*J(:,j+1);
end
GFD=gfd0;

%% ----  2nd
mu=0;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-.6;

Coe=JacobiApprCoe(@gg,0,0,N);
Lx=(max(log(x))-min(log(x)));
X=2*(log(x)-min(log(x)))/Lx-1;
J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
gfd0=zeros(size(x'));
for j=0:N
    gfd0=gfd0+Coe(j+1)*J(:,j+1);
end
GFD=[GFD,gfd0];

%% ---- 3rd
mu=1;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-1;

Coe=JacobiApprCoe(@gg,0,0,N);
Lx=(max(log(x))-min(log(x)));
X=2*(log(x)-min(log(x)))/Lx-1;
J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
gfd0=zeros(size(x'));
for j=0:N
    gfd0=gfd0+Coe(j+1)*J(:,j+1);
end
GFD=[GFD,gfd0];

%% ---- 4th
mu=1;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-.6;

Coe=JacobiApprCoe(@gg,0,0,N);
Lx=(max(log(x))-min(log(x)));
X=2*(log(x)-min(log(x)))/Lx-1;
J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
gfd0=zeros(size(x'));
for j=0:N
    gfd0=gfd0+Coe(j+1)*J(:,j+1);
end
GFD=[GFD,gfd0];
%
% plot(x,GFD,'LineWidth',2,'MarkerSize',10)
% set(gca,'FontSize',12);
% grid on
% xlabel('x')
% ylabel('Hadamard fractional integral')
% legend('w(x)=1,\alpha=1','w(x)=1,\alpha=0.6','w(x)=x,\alpha=1','w(x)=x,\alpha=0.6');

%% --------------------------------------------

NN=8:22;

%%  ---- 1st
mu=0;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-1;
err=[];
GFD(:,1)=sinint(pi*x)-sinint(pi);
for N=NN
    Coe=JacobiApprCoe(@gg,0,0,N);
    Lx=(max(log(x))-min(log(x)));
    X=2*(log(x)-min(log(x)))/Lx-1;
    J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
    gfd0=zeros(size(x'));
    for j=0:N
        gfd0=gfd0+Coe(j+1)*J(:,j+1);
    end
    err=[err;max(abs(gfd0-GFD(:,1)))];
end
% Err=[Err;err];
semilogy(NN,err,'b-*','LineWidth',1.5,'MarkerSize',8)
hold on

%% ----  2nd
mu=0;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-.6;
err=[];
for N=NN
    Coe=JacobiApprCoe(@gg,0,0,N);
    Lx=(max(log(x))-min(log(x)));
    X=2*(log(x)-min(log(x)))/Lx-1;
    J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
    gfd0=zeros(size(x'));
    for j=0:N
        gfd0=gfd0+Coe(j+1)*J(:,j+1);
    end
    err=[err;max(abs(gfd0-GFD(:,2)))];
end
% Err=[Err;err];
semilogy(NN,err,'r-o','LineWidth',1.5,'MarkerSize',8)
hold on

%% ---- 3rd
mu=1;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-1;
err=[];
for N=NN
    Coe=JacobiApprCoe(@gg,0,0,N);
    Lx=(max(log(x))-min(log(x)));
    X=2*(log(x)-min(log(x)))/Lx-1;
    J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
    gfd0=zeros(size(x'));
    for j=0:N
        gfd0=gfd0+Coe(j+1)*J(:,j+1);
    end
    err=[err;max(abs(gfd0-GFD(:,3)))];
end
% Err=[Err;err];
semilogy(NN,err,'g-s','LineWidth',1.5,'MarkerSize',8)
hold on

%% ---- 4th
mu=1;ww=inline(strcat('x.^',num2str(mu)),'x');
alp=-.6;
err=[];
for N=NN
    Coe=JacobiApprCoe(@gg,0,0,N);
    Lx=(max(log(x))-min(log(x)));
    X=2*(log(x)-min(log(x)))/Lx-1;
    J=(Lx/2)^(-alp)*J0_Frac(X,0,0,N,alp,'RL','L','A');
    gfd0=zeros(size(x'));
    for j=0:N
        gfd0=gfd0+Coe(j+1)*J(:,j+1);
    end
    err=[err;max(abs(gfd0-GFD(:,4)))];
end
% Err=[Err;err];

%% ------- plot error
semilogy(NN,err,'c-^','LineWidth',1.5,'MarkerSize',8)
set(gca,'FontSize',12);
xlabel('x')
ylabel('maximum error')
legend('w(x)=1,\alpha=1','w(x)=1,\alpha=0.6','w(x)=x,\alpha=1','w(x)=x,\alpha=0.6');





    function g=gg(x)
        x=(x/2+1/2)*(log(2)-log(1))+log(1);
        g=ww(invz(x)).*ff(invz(x));
    end
end