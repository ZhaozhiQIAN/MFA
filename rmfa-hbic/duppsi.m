function Psi=duppsi(mix, tS, j, alg, up_ops)
%DUPPSI update \Psi by direct compuation described in [1].
%
%Input: mix, tS: =\bS in ECME; =\btS in CM, up_ops
%ops=1(default): only the last d-i columns of C_i is updated.
%   =0         : all columns of C_i is updated.
%Output: mix with new sfa.psi.


%   reference: [1] ML Estimation for Factor Analysis: EM or non-EM? by Jianhua Zhao,
%   Philip L. H. Yu and Qibao Jiang. Department of Statistics and Actuarial
%   Science, The University of Hong Kong

%	Copyright (c) Jianhua Zhao (2007)
if nargin<5
    up_ops=1;
end
xdim=mix.nin;
etat=mix.eta./mix.psi(j, :)-1;
switch alg
    case {'ECM1', 'ECM2', 'ECM3'}
        C=eye(xdim)+mix.U{j}.*repmat(1./mix.lambda{j}-1, xdim, 1)*mix.U{j}';
    case {'ECME3-1', 'ECME3-2', 'ECME3-3'}
        psisinv=mix.psi.^-.5;
        psiinvA=repmat(psisinv', 1, mix.subdim).*mix.A;
        tS=repmat(psisinv', 1, xdim).*tS.*repmat(psisinv, xdim, 1);
        C=eye(xdim)-psiinvA*inv(eye(mix.subdim)+psiinvA'*psiinvA)*psiinvA';
    otherwise
        error([' Unknown algorithms']);
end
for i=1:xdim
    a=C(i,i);
    Ccol=C(:, i);
    omega=a^-2*(Ccol'*tS*Ccol-a);
    omega=max(omega, etat(i));
    if omega>etat(i) || mix.psi(j, i)~=mix.eta % if old \psi=eta and \omega_i=etat_i, no need to update \psi and C
        r=omega/(1+omega*a);
        if i<xdim
            if up_ops C(:, i+1:xdim)=C(:, i+1:xdim)-r*Ccol*Ccol(i+1:xdim)';else C=C-r*C(:, i)*C(:, i)';end
        end
        mix.psi(j, i)=max((omega+1)*mix.psi(j, i), mix.eta);
    end
end
Psi=mix.psi(j, :);

