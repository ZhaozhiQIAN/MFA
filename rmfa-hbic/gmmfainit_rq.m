% ranking with different q: A formated as in different q (eff_dim * X dim)
function mix = gmmfainit_rq(x, ncentres, subdim, options, start, alg)

mix.ncentres=ncentres;
if length(subdim)==1
    mix.subdim=repmat(subdim, 1, mix.ncentres);
    mix.effdim=mix.subdim;
else
    mix.subdim=subdim;
    mix.effdim=subdim;
end
[ndata, xdim]= size(x);
mix.nin=xdim;
mix.eta=5e-3;

% % use option 1
% [centres_old, e_old, post_old] = kmeans(ncentres,x, options);
% for k=1:2
%      [centres_new, e_new, post_new] = kmeans(ncentres,x, options);
%      if e_new<e_old
%            centres_old=centres_new;
%            post_old=post_new;
%            e_old=e_new;
%      end
% end
mix.centres = mfakmeans(x, ncentres, options);
        
mix.priors=1/mix.ncentres*ones(1, mix.ncentres); 
post=repmat(1/mix.ncentres, ndata, mix.ncentres);
new_pr=sum(post, 1);

%set prior
n2 = dist2(x, mix.centres);
mix.pripara=0.1;%0.1;
sig=sum(sum(post.*n2,1))/(ndata*mix.nin);
mix.pri=mix.pripara*mix.nin*sig;
for j=1:mix.ncentres
    mix.U{j} = eye(mix.nin,mix.subdim(j));
    mix.lambda{j} = sig*ones(1,mix.subdim(j));
    mix.psi(j, :) = repmat(mix.lambda{j}(1)/10,1,xdim);
    mix.A{j}=sqrt(mix.lambda{j}-1)'*sqrt(mix.psi(j, :)).*mix.U{j}';
end

mix.nwts = 1 + mix.nin*(mix.effdim+2) - mix.effdim.*(mix.effdim-1)/2 ;






