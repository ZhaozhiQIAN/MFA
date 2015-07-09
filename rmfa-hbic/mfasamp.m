function x=mfasamp(mfa, ndata) 
%sample ndata*ncentres points from MFA model

%ndata number of points in sub-model

% Copyright (c) Jianhua Zhao (2006)
% Dept. of Statistics & Actuarial Science, The University of Hong Kong

for i = 1:mfa.ncentres
    ymean = zeros(mfa.subdim(i), 1);
    ysigma = eye(mfa.subdim(i));
    y{i} = mvnrnd(ymean, ysigma, ndata);
end

for j=1: mfa.ncentres
    start=(j-1)*ndata;
    noisdata=mvnrnd(mfa.centres(:, j), diag(mfa.psi(j, :)), ndata);
    x(:, start+1: start+ndata)=(y{j}*mfa.A{j}'+noisdata)';
end
p = randperm(mfa.ncentres*ndata);
x = x(:,p);
