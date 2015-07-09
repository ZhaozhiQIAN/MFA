% ranking with different q: A formated as in different q (eff_dim * X dim)
function mix = gmmfainit_rq(x, ncentres, subdim, options, ~, ~)

subdim = max(subdim);

mix = gmmfainit(x, ncentres, subdim, options, 'pca', 'ECM2');

for i = 1:ncentres
    mix.A{i} = mix.A{i}';
end
mix.priors = mix.priors';
mix.centres = mix.centres';
% updated
mix.nwts = 1 + (mix.nin-1)*(mix.effdim+2) - mix.effdim.*(mix.effdim-1)/2 -1;
mix.pripara=0.1;
mix.covar_type = 'mfa';
