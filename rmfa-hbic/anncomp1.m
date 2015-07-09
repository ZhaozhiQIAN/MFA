function [mix, post]=anncomp1(x, mix, post)
%merge the component with the least prob with another one with least kullback distance
[ndata, xdim] = size(x);
[minpi, idx]=min(mix.priors);D=kulldist(mix,idx);[minD,ind]=min(D);
% merge ind and idx
cand=[ind idx];
%elipsmix(mix,2)
% D
% ndata*mix.priors
% mix.pri
% update prior prob and center
mix.priors(ind)=sum(mix.priors(cand));
mix.centres(ind,:)=mix.priors(cand)*mix.centres(cand,:)/mix.priors(ind);
post(ind,:)=sum(post(cand,:),1);
mix=delcomq(mix, idx); post(idx,:)=[];
% update parameter prior
n2 = dist2(x, mix.centres);
mix.pri=mix.pripara*mix.nin*sum(sum(post'.*n2,1))/(ndata*mix.nin);

function D=kulldist(mix,idx)
if mix.effdim(idx)>0
    psiinv=1./mix.psi(idx, :);
    psiinvA=mix.A{idx}.*repmat(psiinv, mix.effdim(idx), 1);
    Minv=inv(eye(mix.effdim(idx))+mix.A{idx}*psiinvA');
    d0=sum(log(mix.psi(idx, :)))-log(det(Minv));
else
    d0=sum(log(mix.psi(idx, :)));
end
for j=1:mix.ncentres
    d=zeros(1,6);
    if j~=idx
        if mix.effdim(j)>0
            psiinv=1./mix.psi(j, :);
            psiinvA=mix.A{j}.*repmat(psiinv, mix.effdim(j), 1);
            Minv=inv(eye(mix.effdim(j))+mix.A{j}*psiinvA');
             if mix.effdim(idx)>0
            d(1)=sum(sum(mix.A{idx}.*repmat(psiinv, mix.effdim(idx), 1).*mix.A{idx}));
            psiinvAA0=psiinvA*mix.A{idx}';
            d(2)=sum(sum(Minv*psiinvAA0.*psiinvAA0));
             end
            d(3)=sum(psiinv.*mix.psi(idx,:));
            psiinvApsi0=psiinvA.*repmat(mix.psi(idx,:), mix.effdim(j), 1);
            d(4)=sum(sum(Minv*psiinvApsi0.*psiinvA));
            mu=mix.centres(idx, :)-mix.centres(j, :);
            mupsi=mu.*psiinv;muA=psiinvA*mu';muAM=Minv*muA;
            d(5)=mupsi*mu'-muAM'*muA;
            d(6)=sum(log(mix.psi(j, :)))-log(det(Minv))-d0;
            D(j)=sum(d)-mix.nin;
        else
            psiinv=1./mix.psi(j, :);      
            if mix.effdim(idx)>0
            d(1)=sum(sum(mix.A{idx}.*repmat(psiinv, mix.effdim(idx), 1).*mix.A{idx}));         
            end
            d(2)=sum(psiinv.*mix.psi(idx,:));           
            mu=mix.centres(idx, :)-mix.centres(j, :);
            mupsi=mu.*psiinv;
            d(3)=mupsi*mu';
            d(4)=sum(log(mix.psi(j, :)))-d0;
            D(j)=sum(d)-mix.nin;           
        end
    else
        D(j)=inf;
    end
end






