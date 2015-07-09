% Return Values
%   mix: mix.priors ( Prob data come from cluster i) and mix.centres (mu) are updated
%
%   S: local covariance matrix for each component (stored in 3D array, dim 3 = j)
%
%   logl: log likelihood
%   post: expected membership variables
% Known Issues
%   When dimension is too big, mix.centres, S, logl might be NA.

function [mix, S, logl, postc]=estep1(r, mix, I, u, ~)
%increase no. of draws if L^(t+1)<L^t
%post: Mxn, x1_t:dxM
ndraw=size(u, 2);
ndata=sum(r.c);[xdim, npat]= size(r.p);
S=I.C0;
sigma=I.C0;J.e=I.e(1:ndraw); J.t=I.t(:, 1:ndraw);
logP=I.p;loga=logP(:, 1); post=logP; postc=logP;
x1_t=I.x1_t; x2_t=I.C0; mu_diff1=x1_t;
% resolve conflict (row/column)
centres = mix.centres';
priors = mix.priors';

logPro=log(priors+realmin);

for j=1:mix.ncentres
    sigma(:, :, j)=mix.A{j}'*mix.A{j}+diag(mix.psi(j, :));
end

for i=1:npat
    for j=1:mix.ncentres
        ur=rhaltvec(u);
        ind=r.p(:, i);
        [tmp,ind2]=sort(ind);
        % calling ghk simulator
        %	OUTPUT:	p = ranking prob.
        %           y.y1 = E(x); 
        %           y.y2 = E(xx'); 
        [y, p]=ghkr(I.C*centres(ind,j), I.C*sigma(ind, ind, j)*I.C',ur, J);
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
new_pr=sum(postc, 2);
priors=new_pr/ndata;
mu_diff1=mu_diff1./repmat(new_pr', xdim, 1);
mu_diff2=mu_diff1;
%mu_diff2(end, :)=zeros(1, mix.ncentres);
centres=mu_diff2+centres;
for j=1:mix.ncentres
    S(:, :, j)=S(:, :, j)/new_pr(j)-mu_diff1(:, j)*mu_diff2(:, j)'-mu_diff2(:, j)*mu_diff1(:, j)'+mu_diff2(:, j)*mu_diff2(:,j)';
end
logl= sum(sum((-log(post+realmin) + logP).*postc,2));
mix.centres = centres';
mix.priors = priors';
% scale and location shift
scl = 1:mix.ncentres;
for j=1:mix.ncentres
    scl(j) = mean(diag(S(:, :, j)));
%     scl = max(max(S(:, :, j)));
    S(:, :, j)=S(:, :, j)/scl(j);
    mix.centres(j,:) = mix.centres(j,:)/sqrt(scl(j));
    mix.centres(j,:) = mix.centres(j,:) - mean(mix.centres(j,:));
    mix.psi(j, :) = mix.psi(j, :)/scl(j);
end

return

