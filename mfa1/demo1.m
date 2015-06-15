%This M-file is used to compare different ML algorithms for fitting MFA.
clear all;
clf;
%% generate sample
dim=7;ndata=1000;dppca_dim=2;
mfa.ncentres=2;      mfa.A=zeros(dim, dppca_dim, mfa.ncentres);
mfa.psi=zeros(mfa.ncentres, dim);    mfa.centres=0.1*[1:dim;dim:-1:1]';
for j=1: mfa.ncentres
    mfa.A(:, :, j)=3*rand(dim, dppca_dim);%x5 for GEM2 
    mfa.psi(j, :)=j*rand(1, dim); %mfa.psi(j, 2:20:dim)=200*rand(1, length(2:20:dim));
end
% mfasamp: sample data from mfa model
mfa.subdim=dppca_dim; x=mfasamp(mfa, ndata);
%% fit the data
clf;
[y, ro]=sort(x, 1, 'ascend');
profile on;
g=group(ro');
r.p=ro(:, g.pi);
r.c=g.c;

options(1)=0;options(3)=1e-4;options(14)=20;
ncentres=mfa.ncentres;    subdim=mfa.subdim; start='pca';

alg={'ECM2' 'EM'};    colr={'b' 'g'};
% fitting options are stored in options
options(1)=0;options(3)=1e-1; 
options(14)=100;options(20)=50;options(21)=200;options(22)=2000;
% initialize parameters inside mix1
mix1 = gmmfainit(ro', ncentres, subdim, options, start, 'ECM2');
for i=1:2
    tic
    if strcmp(alg{i}, 'EM')
        %convert mix used in ECM to be used in EM 
        mix1=cstru(mix1, 'ECM2');
    end
    [mix, options, errlog, ndraw]=gmmfa1(mix1, r, options, alg{i});%,
    T(i)=toc;
    iternum(i)=options(16);
    LogL(i)=errlog(options(16));
    ll{i}=errlog(1:options(16));%/(ndata*ncentres);
    lengthll=length(ll{i});
    h(i)=plot(ll{i}, strcat(colr{i}, '-'));
    hold on;
    ndw{i}=ndraw(1:iternum(i));
end
legend(h,alg);
profile off;
profile viewer;