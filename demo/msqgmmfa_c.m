function [mix, options, errlog,post,logact,kill] = msqgmmfa_c(mix, x, options, method)
%GMMFA	EM algorithms for Gaussian mixture model.
% criterion value is stored in options(8)
%x, x2: dxn; Ax=qxnxM; A:qxd; psiinvA:qxd; loga: MxN; post: MxN;

[xdim, ndata] = size(x);
% Sort out the options
if options(14) niters = options(14);else niters = 100;end
display = options(1);store = 0;
if (nargout > 2) store = 1;	% Store the error values to return them
    errlog = zeros(1, niters);end
test = 0;
if options(3)
    test = 1;
end
% Test: test log likelihood for termination
% L_sig: signal the iteration which satisfied a weaker convergence
eold=-inf; kill=[];
x2=x.^2; xt=x';
n0=15;   oldncen=mix.ncentres;

% Main loop of algorithm
for n = 1:niters
    %  E-step
    [post, logact,logL]= gmmfapost_c(mix, x, x2); %E-step
    switch method
        case 'HBIC'
            % e is HBIC value
            e = logL - sum((mix.nwts-1)/2.*log(max(ndata*mix.priors,1)))...
                - (mix.ncentres-1)/2*log(ndata);
%             logL
%             sum((mix.nwts-1)/2.*log(max(ndata*mix.priors,1)))- (mix.ncentres-1)/2*log(ndata)
            %                         if (oldncen==mix.ncentres & e<eold)
            %                             mix.ncentres
            %                             e-eold
            %                             fprintf('HBIC decreases\n');
            %                         end
            %     new_pr = sum(post, 2);
        case 'MML'
            e = logL - sum((mix.nwts-1)/2.*log(max(ndata*mix.priors,1)))...
                - mix.ncentres/2*log(ndata)+sum(mix.nwts)/2*(log(12)-1);
        case 'ML'
            e=logL;
        otherwise
            error(['Unknown criterion ', method]);
    end
    % CM-step 1: update \pi_j (marginal probability of cluster)
    mix.priors = (sum(post, 2)/ ndata)';
    oldncen=mix.ncentres;
    % Annihilate component
    switch method
        case {'HBIC','MML'}
            while min(ndata*mix.priors)<=n0
                [mix, post]=anncomp1(xt, mix, post);
            end
        case 'ML'
            while min(ndata*mix.priors)<=1
                [mix, post, logact]=anncomp(xt, mix, post, logact, 1);
            end
        otherwise
            error(['Unknown criterion ', method]);
    end
    if mix.ncentres<oldncen
        fprintf('annihilate %d components\n', oldncen-mix.ncentres);
        kill=[kill n];
    end
    % CM-step 1(continued): update \pi_j after annihiliation (stored in
    % new_pr)
    new_pr = sum(post, 2);
    for j=1: mix.ncentres
        %  membership variable R is post
        %    mix.nin = dim of x
        Rx=repmat(post(j, :), mix.nin, 1).*x;%dxN
        %  CM-step 1: update \mu_j
        mix.centres(j, :) = sum(Rx, 2)/new_pr(j); 
        %  calculate Sj
        S=(Rx*x'-mix.centres(j, :)'*mix.centres(j, :)*new_pr(j)+mix.pri*eye(mix.nin))/(new_pr(j)+mix.nin);
        S1=(S+S')/2; %guarantee symetric S; S1=S;
        psisinv=mix.psi(j, :).^-0.5;
        S=psisinv'*psisinv.*S1;
        %  CM-step 1.5: select q (q is stored in mix.effdim)
        switch method
            case 'HBIC'
                [tempU, templambda, mix.effdim(j)] = msfa(S, new_pr(j), new_pr(j));
            case 'MML'
                [tempU, templambda, mix.effdim(j)] = msfa(S, new_pr(j), new_pr(j)/12*exp(1));
            case 'ML'
                [templambda,tempU] = eigdec(S, mix.subdim(j));
                templambda=templambda(templambda>1)';
                mix.effdim(j)=length(templambda);
            otherwise
                error(['Unknown criterion ', method]);
        end
        if mix.effdim(j)>0
            mix.lambda{j}=templambda;           
            mix.U{j}=tempU(:, 1:mix.effdim(j));
            %  CM-step 2: Update A
            mix.A{j}=sqrt(mix.lambda{j}-1)'*sqrt(mix.psi(j, :)).*mix.U{j}';
            %  CM-step 3: Update \Psi
            mix.psi(j, :)=uppsi(mix, S, j);         
        else
            mix.lambda{j}=[];
            mix.U{j}=[];
            mix.A{j}=[];
            mix.psi(j, :)=max(diag(S1)',mix.eta);
        end
    end
    mix.nwts = 1 + mix.nin*(mix.effdim+2) - mix.effdim.*(mix.effdim-1)/2 ;
    options(16)=n;
    
    if (display || store || test)
        if store  errlog(n) = e;  end
        if (n>1 && display)   %fprintf(1, 'Cycle %4d  logL %11.6f, relative increment %e\n', n, e, 1-abs(e/eold));
            if (n > 1 && e<eold)  fprintf('----> LogL decreased in iteration %4d\n', n);
                fprintf(1, 'Cycle %4d  logL %11.6f, relative increment %e\n', n, e, 1-abs(e/eold));
            end
        end
        if test
            if (n > 1 && abs(1 - eold/e) < options(3))
                options(8) = e; return;
            else  eold = e; end
        end
    end
end
options(8) = errlog(n);
fprintf(1, 'Cycle %4d  logL %11.6f\n', n, options(8));
if display disp('Warning: Maximum number of iterations has been exceeded');end
