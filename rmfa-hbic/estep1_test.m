%% setting up test cases
cd D:\MFA\rmfa-hbic
% set "dppca_dim=[2,2];" in sample_ranking.m
sample_ranking

options(1)=0;options(3)=1e-4;options(14)=20;
ncentres=9;    subdim=2; start='pca';

% fitting options are stored in options
options(1)=0;options(3)=1e-1; 
options(14)=100;options(20)=50;options(21)=200;options(22)=2000;
% initialize parameters inside mix1
cd D:\MFA\mfa1
mix1 = gmmfainit(ro', ncentres, subdim, options, start, 'ECM2');

ndraw_max=options(22);
ndraw(1:2)=options(20);
ndrawi=options(21);
[xdim, npat]= size(r.p);
ndata=sum(r.c);
C0=zeros(xdim);I.C=C0;
for j=1:xdim-1
    I.C(j, j: j+1)=[1 -1];
end
I.C(xdim, xdim)=1;I.Ci=inv(I.C);I.C0=zeros(xdim, xdim, mix1.ncentres);
I.p=zeros(mix1.ncentres,npat);I.x1_t=zeros(xdim,mix1.ncentres);
I.e=ones(1, ndraw_max);I.t=zeros(xdim-1, ndraw_max);
u=haltvec(xdim-1, ndraw_max);

tol = 0.05;
%% test one: original case
[mix_org, S_org, e_org]=estep1(r, mix1, I, u(:, 1:ndraw_max), 'ECM2');

cd D:\MFA\rmfa-hbic
for i = 1:ncentres
    mix1.A{i} = mix1.A{i}';
end

mix1.priors = mix1.priors';
mix1.centres = mix1.centres';
[mix_test, S_test, logl_test, post]=estep1(r, mix1, I, u(:, 1:ndraw_max));

assert(isempty(find(abs(mix_test.priors-mix_org.priors')>tol, 1)));
assert(isempty(find(abs(S_org-S_test)>tol, 1)));
assert(isempty(find(abs(mix_test.centres-mix_org.centres')>tol, 1)));
assert(~isnan(mix_test.centres))
disp('test1 finish')
%% test two: randomly assign initial A
A_rand1 = rand(15,10);
A_rand2 = rand(15,10);
mix1.A{1} = A_rand1;
mix1.A{2} = A_rand2;
mix1.priors = mix1.priors';
mix1.centres = mix1.centres';
cd D:\MFA\mfa1
[mix_org, S_org, e_org]=estep1(r, mix1, I, u(:, 1:ndraw_max), 'ECM2');
cd D:\MFA\rmfa-hbic
mix1.A{1} = mix1.A{1}';
mix1.A{2} = mix1.A{2}';
mix1.priors = mix1.priors';
mix1.centres = mix1.centres';
[mix_test, S_test, logl_test, post]=estep1(r, mix1, I, u(:, 1:ndraw_max));

assert(isempty(find(abs(mix_test.priors-mix_org.priors')>tol, 1)));
assert(isempty(find(abs(S_org-S_test)>tol, 1)));
assert(isempty(find(abs(mix_test.centres-mix_org.centres')>tol, 1)));
assert(~isnan(mix_test.centres))
disp('test2 finish')

%% test three: randomly assign initial Psi
Psi_rand = rand(ncentres,15);
mix1.psi = Psi_rand;
mix1.A{1} = mix1.A{1}';
mix1.A{2} = mix1.A{2}';
mix1.priors = mix1.priors';
mix1.centres = mix1.centres';
cd D:\MFA\mfa1
[mix_org, S_org, e_org]=estep1(r, mix1, I, u(:, 1:ndraw_max), 'ECM2');
cd D:\MFA\rmfa-hbic
mix1.A{1} = mix1.A{1}';
mix1.A{2} = mix1.A{2}';
mix1.priors = mix1.priors';
mix1.centres = mix1.centres';
[mix_test, S_test, logl_test, post]=estep1(r, mix1, I, u(:, 1:ndraw_max));

assert(isempty(find(abs(mix_test.priors-mix_org.priors')>tol, 1)));
assert(isempty(find(abs(S_org-S_test)>tol, 1)));
assert(isempty(find(abs(mix_test.centres-mix_org.centres')>tol, 1)));
assert(~isnan(mix_test.centres))

disp('Passed all three test cases')













