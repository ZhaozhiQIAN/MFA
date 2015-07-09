function mix = delcomq(mix, l)
mix.ncentres=mix.ncentres-length(l);
mix.priors(l)=[];
mix.centres(l, :)=[];
mix.lambda(l)=[];
mix.psi(l,:)=[];
mix.U(l)=[];
mix.A(l)=[];
if strcmp(mix.covar_type,'tmfa')
    mix.v(l)=[];
end
mix.effdim(l)=[];
mix.subdim(l)=[];
mix.nwts(l)=[];
