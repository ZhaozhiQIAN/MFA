% target: 
%   Show if using eigenvalues (equation 32) results in large q.
% methodology: 
%   Perform a grid search on all combinations of q. For each configuration
%   calculate log-likelihood using numeric method. Select model using HBIC
%   criterion (equation 29). Due to large amount of computation, two mixing
%   components are used in the experiment. Repeat data generation and
%   estimation for 6 times.
% observation:
%   EM converges within 10 steps (relative change of HBIC is less than
%   1e-3). The selected q values are close to true values.
% reproducible example:

rng(int32(37));
% parameters: the same as demo1.m
dim=10;
ndata=1000;
dppca_dim=[2,2];
mfa.ncentres=2;      
mfa.psi=zeros(mfa.ncentres, dim);    
mfa.centres=0.1*[1:dim;dim:-1:1]';
mfa.A{1} = 3*rand(dim,dppca_dim(1));
mfa.A{2} = 3*rand(dim,dppca_dim(2));
for j=1: mfa.ncentres
    mfa.psi(j, :)=j*rand(1, dim);
end

for run = 1:6
    % mfasamp: sample data from mfa model
    mfa.subdim=dppca_dim; 
    x=mfasamp(mfa, ndata);
    % generate rankings
    [y, ro]=sort(x, 1, 'ascend');
    g=group(ro');
    % r stores rankings: r.p is ranking format (column wise), r.c is count vector 
    r.p=ro(:, g.pi);
    r.c=g.c;
    %% fitting model
    % set options
    options(1)=0; % display
    options(3)=1e-3; 
    options(14)=30;options(20)=50;options(21)=200;options(22)=2000;
    % loop over ncentres
    cri_hbic_best = -Inf;
    %cri_mml_best = -Inf;
    for ncentres=2:2
        subdim=repmat(6,1,ncentres); 
        start='pca';
        options(1)=0;
        mix1 = gmmfainit_rq(ro', ncentres, subdim, options, start, 'ECM2');
        options(1)=1;
        [mix_hbic{ncentres}, options_hbic{ncentres}, errlog_hbic{ncentres}, ndraw_hbic{ncentres}]=rank_em_q_main(mix1, r, options, 'HBIC');
        %[mix_mml{ncentres}, options_mml{ncentres}, errlog_mml{ncentres}, ndraw_mml{ncentres}]=rank_em_q_main(mix1, r, options, 'MML');
        if options_hbic{ncentres}(8) > cri_hbic_best
            mix_hbic_best{run} = mix_hbic{ncentres};
            cri_hbic_best = options_hbic{ncentres}(8);
        end
        %if options_mml{ncentres}(8) > cri_mml_best
        %    mix_mml_best{run} = mix_mml{ncentres};
        %    cri_mml_best = options_mml{ncentres}(8);
        %end
    end;
end;


for i=1:run
    fprintf('>>>> Run %d <<<<\n',i);
    disp('hbic')
    disp(mix_hbic_best{run}.ncentres)
    for c=1:mix_hbic_best{run}.ncentres
        fprintf('\tcomponent %d: q = %d\n', c, size(mix_hbic_best{run}.A{c},1));
    end
%     disp('mml')
%     disp(mix_mml_best{run}.ncentres)
%     for c=1:mix_mml_best{run}.ncentres
%         fprintf('\tcomponent %d: q = %d\n', c, size(mix_mml_best{run}.A{c},1));
%     end
end




