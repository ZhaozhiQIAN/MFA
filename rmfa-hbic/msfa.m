function [U, lambda, subdim] = msfa(x, ndata1, ndata2, regpara)
%PPCA	Probabilistic Principal Components Analysis
%
%	Description
%	 [VAR, U, LAMBDA] = PPCA(X, PPCA_DIM) computes the principal
%	component subspace U of dimension PPCA_DIM using a centred covariance
%	matrix X. The variable VAR contains the off-subspace variance (which
%	is assumed to be spherical), while the vector LAMBDA contains the
%	variances of each of the principal components.  This is computed
%	using the eigenvalue and eigenvector  decomposition of X.
%
%	See also
%	EIGDEC, PCA
%

%	Copyright (c) Jianhua Zhao (1996-2001)
% regpara=1e-1;
data_dim = size(x, 2);
% Assumes that x is centred and responsibility weighted
% covariance matrix
[l Utemp] = eigdec(x, data_dim);
% Zero any negative eigenvalues (caused by rounding)
l(l<0) = 0;
qmax=floor(data_dim+0.5*(1-sqrt(1+8*data_dim)));%=15;%
q_temp = min(max(find(l> 1)),qmax);
%q_temp = data_dim-min(find(s2_temp> 1e-2));
if length(q_temp)==0 %length(q_temp) == 0 || q_temp==0
    % All the latent dimensions have disappeared, so we are
    % just left with the noise model    
    lambda = 0;%var;
    U = [];%Utemp(:, 1);
    subdim=0;%ppca_dim=1;
else    
    npara=(data_dim-1)*([1: q_temp]+2)-[1: q_temp].*([1: q_temp]-1)/2 - 1;    
   % ind=[data_dim-1:-1:data_dim-q_temp];
   % BIC=-ndata1*(cumsum(log(l(1: q_temp)')) + ind.*log(s2_temp(ind)'))-npara*log((ndata2));
    BIC=-ndata1*(cumsum(log(l(1: q_temp)')-l(1: q_temp)'+1))-npara*log((ndata2));
    % [y, ppca_dim]=max(BIC);
    %[y, ppca_dim]=max(-ndata1*(cumsum(log(l(1: q_temp)')) + ind.*log(s2_temp(ind)'))-npara*log((ndata2)));
    %fist augument
    tmp=find(BIC>=[BIC(2:end) -1e+10]);
%     if numel(tmp)==0
%         tmp=1;
%     end
    subdim=tmp(1);
    %    [y,ppca_dim]=max(-ndata1/2*(data_dim*log(s2_temp(ind)')+ind(end:-1:1).*(1+log(2*pi))));%-npara/2*(log(ndata2)+1));%
    U = Utemp(:, 1:subdim);
    lambda= l(1:subdim)';
end



