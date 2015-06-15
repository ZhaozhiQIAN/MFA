function loga = gmmfaactiv_c(mix, r, options)
%GMMACTIV Computes the activations of a Gaussian mixture model.
% output: Minv, psiinvA, xA to avoid re-computation in M-step.
%x, x2: dxn; Ax=qxnxM; A:qxd
ndraw=options(15);
ndata = size(x, 2); loga = zeros(mix.ncentres, ndata);  % Preallocate matrix
log_normal = mix.nin*log(2*pi); d2 = zeros(mix.ncentres, ndata);
logZ = zeros(mix.ncentres, 1); cal=0; if nargin<5 cal=1;end

for j = 1:mix.ncentres
    psiinv=1./mix.psi(j, :);
    psiinvA=mix.A{j}.*repmat(psiinv, mix.effdim(j), 1);
    Minv=inv(eye(mix.effdim(j))+mix.A{j}*psiinvA');
    Ax=psiinvA*x;
    mupsi=mix.centres(j, :).*psiinv;
    muA=psiinvA*mix.centres(j, :)';
    muAM=Minv*muA;
    logZ(j)=log_normal+sum(log(mix.psi(j, :)))-log(det(Minv))+mupsi*mix.centres(j, :)'-muAM'*muA;
    d2(j, :)= psiinv*x2 + (-2*mupsi)*x - sum(Minv*Ax.*Ax, 1)+(2*muAM)'*Ax;
    ind=find(d2(j, :)+repmat(mupsi*mix.centres(j, :)'-muAM'*muA, 1, ndata)<0);
    if length(ind)>0
        d2(j, ind)+mupsi*mix.centres(j, :)'-muAM'*muA
    end
end
loga = -0.5*(d2 + repmat(logZ, 1, ndata));



