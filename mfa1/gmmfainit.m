% for rankings
function mix = gmmfainit(x, ncentres, subdim, options, start, alg)
%GMMINIT Initialises Gaussian mfa from data

%initialize parameters
[ndata, xdim]= size(x);
mix.subdim=subdim; %all components share a same subspace dimension.
mix.ncentres=ncentres;mix.nin=xdim;mix.eta=1e-2;
mix.effdim=repmat(subdim, 1, mix.ncentres);

mix.centres = mfakmeans(x, ncentres, options);
% cluster_sizes = max(sum(post, 1), 1);  % Make sure that no prior is zero
% switch start
%     case 'pca'
%         mix.priors = cluster_sizes/sum(cluster_sizes); % Normalise priors
%     case 'adhoc'
%         mix.priors=1/mix.ncentres*ones(1, mix.ncentres);
% end
mix.priors=repmat(1/mix.ncentres,mix.ncentres, 1);

for j = 1:mix.ncentres
    % Pick out data points belonging to this centre
   diffs = x - repmat(mix.centres(j, :), ndata, 1);
    switch alg
        case {'ECM1', 'ECM2'}
            switch start
                case 'pca'
                    [tempcovars, tempU, templambda] = ...
                        ppca((diffs'*diffs)/ndata, mix.subdim);
                    if length(templambda) ~= mix.subdim
                        error('Unable to extract enough components');
                    else
                        mix.U{j}=tempU;
                        mix.lambda{j}=templambda;
                        mix.psi(j, :)=repmat(tempcovars, 1, xdim);
                        mix.A{j}=mix.U{j}.*repmat(sqrt(templambda-tempcovars), xdim, 1);
                    end
                case 'adhoc'
                    mix.U{j}=eye(mix.nin, mix.subdim);
                    mix.lambda{j}=ones(1, mix.subdim);
                    mix.A{j}=mix.U{j}';
                    mix.psi(j, :)=10*ones(1, mix.nin);
            end
        case 'EM'
            switch start
                case 'pca'
                    [tempcovars, tempU, templambda] = ...
                        ppca((diffs'*diffs)/ndata, mix.subdim);
                    if length(templambda) ~= mix.subdim
                        error('Unable to extract enough components');
                    else
                        mix.A(:, :, j)=tempU.*repmat(sqrt(templambda-tempcovars), xdim, 1);
                        mix.psi(j, :)=repmat(tempcovars, 1, xdim);
                    end
                case 'adhoc'
                    mix.A(:, :, j)=eye(mix.nin, mix.subdim);
                    mix.psi(j, :)=10*ones(1, mix.nin);
            end
    end
end
mix.centres=mix.centres';
%mix.centres(end, :)=zeros(1, mix.ncentres);
% mix.psi(:, end)=ones(mix.ncentres, 1);
% mix.priors=mix.priors';