
zn=2/3;
wn=1;
ALP=[0.3 0.6 0.9];
NN=8:36;
for p=6:7
    ERR=zeros(3,29);
    for k=1:3
        alpha=ALP(k);
        for m=1:29
            [rx,M,U,exact]=SpectralGFODEEx41(NN(m),alpha,zn,wn,p);
            ERR(k,m)=max(abs(U-exact));
        end
    end
    figure;
    loglog(NN,ERR(1,:),'b-*',NN,ERR(2,:),'g-*',NN,ERR(3,:),'r-*')
    xlim([8,36])
end