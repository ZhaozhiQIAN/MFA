function [mix, S, logl]=estep1(r, mix, I, u, alg)
%increase no. of draws if L^(t+1)<L^t
%post: Mxn, x1_t:dxM
ndraw=size(u, 2);
ndata=sum(r.c);[xdim, npat]= size(r.p);
S=I.C0;sigma=I.C0;J.e=I.e(1:ndraw); J.t=I.t(:, 1:ndraw);
logP=I.p;loga=logP(:, 1); post=logP; postc=logP;
x1_t=I.x1_t; x2_t=I.C0; mu_diff1=x1_t;

logPro=log(mix.priors+realmin);
switch alg
    case {'ECM1','ECM2','ECM3'}
        for j=1:mix.ncentres
            sigma(:, :, j)=mix.A{j}*mix.A{j}'+diag(mix.psi(j, :));
        end
    case 'EM'
        for j=1:mix.ncentres
            sigma(:, :, j)=mix.A(:, :, j)*mix.A(:, :, j)'+diag(mix.psi(j, :));
        end
end
for i=1:npat
    for j=1:mix.ncentres
        ur=rhaltvec(u);
        ind=r.p(:, i);
        [tmp,ind2]=sort(ind);
        [y, p]=ghkr(I.C*mix.centres(ind,j), I.C*sigma(ind, ind, j)*I.C',ur, J);
        x1_t(:, j)=I.Ci*y.y1;
        x2_t(:, :, j)=I.Ci*y.y2*I.Ci';
        loga(j)=log(p);
    end
    loga=loga+logPro;
    logP(:, i)=loga;
    loga=exp(loga-max(loga));
    post(:, i)=loga/sum(loga);
    postc(:, i)=post(:, i)*r.c(i);
    for j=1:mix.ncentres
        x1_t(:, j)=postc(j, i)*x1_t(:, j);
        mu_diff1(:, j)=mu_diff1(:, j)+x1_t(ind2, j);
        x2_t(:, :, j)=postc(j, i)*x2_t(:, :, j);
        S(:, :, j)=S(:, :, j)+x2_t(ind2, ind2, j);
    end
end
new_pr=sum(postc, 2);mix.priors=new_pr/ndata;
mu_diff1=mu_diff1./repmat(new_pr', xdim, 1);
mu_diff2=mu_diff1;
%mu_diff2(end, :)=zeros(1, mix.ncentres);
mix.centres=mu_diff2+mix.centres;
for j=1:mix.ncentres
    S(:, :, j)=S(:, :, j)/new_pr(j)-mu_diff1(:, j)*mu_diff2(:, j)'-mu_diff2(:, j)*mu_diff1(:, j)'+mu_diff2(:, j)*mu_diff2(:,j)';
end
logl= sum(sum((-log(post+realmin) + logP).*postc,2));
return

