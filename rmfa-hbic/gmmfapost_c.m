function [post, loga, LogL, Minv, psiinvA, Ax, psix2_Ax2] = gmmfapost_c(mix, x, x2, Minv, psiinvA, Ax, psix2_Ax2)
% output: Minv, psiinvA, xA to avoid re-computation in M-step.
%x, x2: dxn; Ax=qxnxM; A:qxd; psiinvA:qxd; loga: MxN; post: MxN;

ndata = size(x, 2);

loga = gmmfaactiv_c(mix, x, x2, 'ECM2');

% post = repmat(mix.priors', 1, ndata).*exp(loga);
% s = sum(post, 1);
% if nargout>1   LogL=sum(log(s+realmin)); end
% % Set any zeros to one before dividing
% s = s + (s==0);
% post = post./repmat(s, mix.ncentres, 1);

% loga_pi = repmat(log(mix.priors+realmin), ndata, 1)+loga;
% post=zeros(ndata, mix.ncentres);
% for i=1:mix.ncentres
%     post(:,i)=1./sum(exp(loga_pi-loga_pi(:,i)*ones(1,mix.ncentres)),2);
% end

LogPxc = loga + repmat(log(mix.priors'+realmin), 1, ndata);%MxN

post    = LogPxc - repmat(max(LogPxc,[],1),mix.ncentres,1);
post    = exp(post)+realmin;
post    = post ./ repmat(sum(post,1),mix.ncentres,1);

if nargout>2 LogL = sum(sum(post.*(-log(post+realmin) + LogPxc),1));end