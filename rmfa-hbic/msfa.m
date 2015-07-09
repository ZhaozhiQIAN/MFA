function q = msfa(e_total, ndata1, ndata2)
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

[q1_dim, q2_dim] = size(e_total);
hbic = zeros(q1_dim, q2_dim);
for q1 = 1:q1_dim
    for q2 = 1:q2_dim
        npara1=(10-1)*(q1+2)-q1*(q1-1)/2 - 1; 
        npara2=(10-1)*(q2+2)-q2*(q2-1)/2 - 1; 
        hbic(q1,q2) = e_total(q1, q2)-npara1*log(ndata1)/2-npara2*log(ndata2)/2;
    end
end

[v,ind]=max(hbic(:));
[q1,q2] = ind2sub(size(hbic),ind);

q = [q1, q2];
