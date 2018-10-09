function [rx,H,Ga]=SolveHadamardIntegral(g,a,b,mu,alp,N)

[rx,r,D]=GFracMat(log(a),log(b),N,alp,@w,@invz);

Ga=a^mu*g(a)/gamma(1-alp)*rx.^(-mu).*(log(rx/a)).^(-alp);
G=g(rx);
H=Ga+D*G;

    function res=w(x)
        res=x.^mu;
    end
    function res=invz(x)
        res=exp(x);
    end
end