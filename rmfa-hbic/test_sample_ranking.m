% unit test for sample_ranking.m
cd D:\MFA\rmfa-hbic
% set "dppca_dim=[2,2];" in sample_ranking.m
sample_ranking
% for backward compatibility
dppca_dim = 2;
mfa.subdim=dppca_dim; 
cd D:\MFA\mfa1

% set options
options(1)=0;
options(3)=1e-4;
options(14)=20;
ncentres=mfa.ncentres;
subdim=mfa.subdim; start='pca';
alg={'ECM2' 'EM'};
options(1)=0;options(3)=1e-1; 
options(14)=100;options(20)=50;options(21)=200;options(22)=2000;
mix1 = gmmfainit(ro', ncentres, subdim, options, start, 'ECM2');

for i=1
    tic
    if strcmp(alg{i}, 'EM')
        %convert mix used in ECM to be used in EM 
        mix1=cstru(mix1, 'ECM2');
    end
    [mix, options, errlog, ndraw]=gmmfa1(mix1, r, options, alg{i});
    T(i)=toc;
    iternum(i)=options(16);
    LogL(i)=errlog(options(16));
    ll{i}=errlog(1:options(16));
    lengthll=length(ll{i});
    hold on;
    ndw{i}=ndraw(1:iternum(i));
end