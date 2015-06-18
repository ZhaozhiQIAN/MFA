function [post, LogL] = gmmfapost_c(mix, x, x2)
%x, x2: dxn; Ax=qxnxM; A:qxd; psiinvA:qxd; loga: MxN; post: MxN;

ndata = size(x, 2);
loga = gmmfaactiv_c(mix, x, x2);
LogPxc = loga + repmat(log(mix.priors'+realmin), 1, ndata);%MxN
post    = LogPxc - repmat(max(LogPxc,[],1),mix.ncentres,1);
post    = exp(post)+realmin;
post    = post ./ repmat(sum(post,1),mix.ncentres,1);
if nargout>1 LogL = sum(sum(post.*(-log(post+realmin) + LogPxc),1));end