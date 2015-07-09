% unit test for sample_ranking.m
cd D:\MFA\rmfa-hbic
% set "dppca_dim=[2,2];" in sample_ranking.m
sample_ranking

% set options
ncentres=mfa.ncentres;
subdim=mfa.subdim; start='pca';
alg={'ECM2' 'EM'};
options(1)=1; % display
options(3)=1e-2; 
options(14)=100;options(20)=50;options(21)=200;options(22)=2000;
mix1 = gmmfainit_rq(ro', ncentres, subdim, options, start, 'ECM2');
mix1.centres = 3 + 0.1*[1:dim;dim:-1:1];
[mix_hbic, options_hbic, errlog_hbic, ndraw_hbic]=rank_em_q_main(mix1, r, options, 'HBIC');
[mix_ml, options_ml, errlog_ml, ndraw_ml]=rank_em_q_main(mix1, r, options, 'ML');

mfa.A{1}*mfa.A{1}'+diag(mfa.psi(1,:));

% for backward compatibility
% dppca_dim = 2;
% subdim=dppca_dim; 
% cd D:\MFA\mfa1
% mix1 = gmmfainit(ro', ncentres, subdim, options, start, 'ECM2');
% mix1.A{1} = mix1.A{1}';
% mix1.A{2} = mix1.A{2}';
% mix1.priors = mix1.priors';
% mix1.centres = mix1.centres';
% 
% [mix_org, options_org, errlog_org, ndraw_org]=gmmfa1(mix1, r, options, 'ECM2');
% cd D:\MFA\rmfa-hbic
