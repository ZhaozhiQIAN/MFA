function mix = gmmfainit1_c(x, ncentres, subdim,covar_type,options)
%gmmfainit1_c Initialises a-mfa

%initialize parameters
[ndata, xdim] = size(x);
mix.ncentres=ncentres; mix.nin=xdim; mix.eta=5e-3; %mix.eta=1e-2;%
mix.covar_type=covar_type;
if length(subdim)==1
    mix.subdim=repmat(subdim, 1, mix.ncentres);
    mix.effdim=mix.subdim;
else
    mix.subdim=subdim;
    mix.effdim=subdim;
end
switch options(5)
    case 2
        perm = randperm(ndata);
        randindex = perm(1:mix.ncentres);
        mix.centres = x(randindex, :);
        %        post=repmat(1/mix.ncentres, ndata, mix.ncentres);
    case 1
        %         opts = statset('MaxIter',20);
        %         [IDX,mix.centres]= kmeans(x, mix.ncentres, 'Replicates',3,'EmptyAction', 'singleton', 'Options',opts);
        [centres_old, e_old, post_old] = kmeans(ncentres,x, options);
        for k=1:2
            [centres_new, e_new, post_new] = kmeans(ncentres,x, options);
            if e_new<e_old
                centres_old=centres_new;
                post_old=post_new;
                e_old=e_new;
            end
        end
        mix.centres=centres_old; post=post_old;  options(8)=e_old;
        %  mix.centres = mfakmeans(x, ncentres, options);
    case 0
        perm = randperm(ndata);
        UCL=min(floor(beta2*ndata), 500);
        perm=perm(1: UCL);
        subx=x(perm, :);
        d2=dist2(subx, subx);
        [Y, col] = max(d2, [], 1);
        [z, row] = max(Y);
        centres=[subx(row, :); subx(col(row), :)];
        [mix.centres, e, post] = kmeans(centres, x, options);
end

mix.priors=1/mix.ncentres*ones(1, mix.ncentres); post=repmat(1/mix.ncentres, ndata, mix.ncentres);
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

switch mix.covar_type
    case 'mfa'
        mix.nwts = 1 + mix.nin*(mix.effdim+2) - mix.effdim.*(mix.effdim-1)/2 ;
    case 'tmfa'
        mix.v=4*one(1,mix.ncentres);
        mix.nwts = 2 + mix.nin*(mix.effdim+2) - mix.effdim.*(mix.effdim-1)/2 ;
end




