Alp=[0.3 0.6 0.9];
cc='brg';
for j=1:3
    alp=Alp(j);
    [rx,H]=SolveHadamardIntegral(@gfun,1,10,1/3,alp,120);
    plot(rx,H,cc(j),'LineWidth',1.5,'MarkerSize',8)
    hold on
end
legend('\alpha=0.3','\alpha=0.6','\alpha=0.9');
xlabel('X')
ylabel('f(X)')
set(gca,'FontSize',12);

