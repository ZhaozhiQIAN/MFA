function x=mfasamp(mfa, ndata) 
%sample ndata*ncentres points from MFA model

%ndata number of points in sub-model

% Copyright (c) Jianhua Zhao (2006)
% Dept. of Statistics & Actuarial Science, The University of Hong Kong


ymean=zeros(mfa.subdim, 1);
ysigma=eye(mfa.subdim);
y=mvnrnd(ymean, ysigma, ndata*mfa.ncentres);
for j=1: mfa.ncentres
    start=(j-1)*ndata;
    noisdata=mvnrnd(mfa.centres(:, j), diag(mfa.psi(j, :)), ndata);
    x(:, start+1: start+ndata)=(y(start+1: start+ndata, :)*mfa.A(:, :, j)'+noisdata)';
end

