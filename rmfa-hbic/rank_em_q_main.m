function [mix, options, errlog,postc,logact,kill] = rank_em_q_main(mix, r, options, method)
%GMMFA	EM algorithms for Gaussian mixture model.
% criterion value is stored in options(8)
%x, x2: dxn; Ax=qxnxM; A:qxd; psiinvA:qxd; loga: MxN; post: MxN;

% Sort out the options
if options(14) niters = options(14);else niters = 100;end
display = options(1);store = 0;
if (nargout > 2) store = 1;	% Store the error values to return them
    errlog = zeros(1, niters);end
test = 0;
if options(3)
    test = 1;
end
ndraw_max=options(22);
ndraw(1:2)=options(20);
ndrawi=options(21);
[xdim, npat]= size(r.p);
ndata=sum(r.c);

% Initialize variables
eold=-inf; 
kill=[];
rt=r.p';
n0=15;   
oldncen=mix.ncentres;
C0=zeros(xdim);
I.C=C0;
for j=1:xdim-1
    I.C(j, j: j+1)=[1 -1];
end
I.C(xdim, xdim)=1;
I.Ci=inv(I.C);
I.C0=zeros(xdim, xdim, mix.ncentres);
I.p=zeros(mix.ncentres,npat);
I.x1_t=zeros(xdim,mix.ncentres);
I.e=ones(1, ndraw_max);
I.t=zeros(xdim-1, ndraw_max);

u=haltvec(xdim-1, ndraw_max);

disp('Start')

% Main loop of algorithm
for n = 1:niters
    if mod(n,10)==0
         fprintf(1, 'Cycle %4d  e %11.6f\n', n, e);
    end
    % calling function estep1
    % update \pi_j (marginal probability of cluster) [mix.priors]
    % update mu [mix.centres], [S]
    % return likelihood [e]
    % update membership [post] for unique rankings; [postc] add number of
    % observations
    last_mix = mix;
    first_time_flag = true;
    time_ind = 0;
    while first_time_flag || isnan(sum(S(:)))
        first_time_flag = false;
        time_ind = time_ind + 1;
        [mix, S, e, postc]=estep1(r, mix, I, u(:, 1:ndraw(n)));
    end 
    if time_ind > 2
        disp(time_ind);
    end
    assert(~isnan(sum(S(:))))
    if (n > 1 && e<eold)
        ndraw(n+1)=min(ndraw(n)+ndrawi, ndraw_max);
    else
        ndraw(n+1)=ndraw(n);
    end

    switch method
        case 'HBIC'
            % e is HBIC value
            e = e - sum((mix.nwts-1)/2.*log(max(ndata*mix.priors,1)))...
                - (mix.ncentres-1)/2*log(ndata);
        case 'MML'
            e = e - sum((mix.nwts-1)/2.*log(max(ndata*mix.priors,1)))...
                - mix.ncentres/2*log(ndata)+sum(mix.nwts)/2*(log(12)-1);
        case 'ML'
            % e = e;
        otherwise
            error(['Unknown criterion ', method]);
    end
    
    oldncen=mix.ncentres;
    % Annihilate component
    switch method
        case {'HBIC','MML'}
            while min(ndata*mix.priors)<=n0
                [mix, postc]=anncomp1(rt, mix, postc);
            end
        case 'ML'
%             while min(ndata*mix.priors)<=1
%                 [mix, postc, logact]=anncomp(rt, mix, postc, logact, 1);
%             end
        otherwise
            error(['Unknown criterion ', method]);
    end
    if mix.ncentres<oldncen
        fprintf('annihilate %d components\n', oldncen-mix.ncentres);
        kill=[kill n];
    end
    
    % CM-step 1(continued): calculate number of data in each component
    component_n = sum(postc, 2);
    assert(sum(component_n)-ndata < 1e-7);
    
    for j=1: mix.ncentres 
        %  CM-step 1.5: select q (q is stored in mix.effdim)
        psisinv=diag(mix.psi(j, :).^-0.5);
        Stil=psisinv*((S(:, :, j)+S(:, :, j)')/2)*psisinv;
        switch method
            case 'HBIC'
                [tempU, templambda, mix.effdim(j)] = msfa(Stil, component_n(j), component_n(j));
            case 'MML'
                [tempU, templambda, mix.effdim(j)] = msfa(Stil, component_n(j), component_n(j)/12*exp(1));
            case 'ML'
                [templambda,tempU] = eigdec(Stil, mix.subdim(j));
                templambda=templambda(templambda>1)';
                mix.effdim(j)=length(templambda);
            otherwise
                error(['Unknown criterion ', method]);
        end
        if mix.effdim(j)>0
            mix.lambda{j}=templambda;           
            mix.U{j}=tempU(:, 1:mix.effdim(j));
            %  CM-step 2: Update A
            if size(mix.A{j},1)~= mix.effdim(j)
                fprintf('Size of A{%d} changed to %d\n',j, mix.effdim(j))
            end
            mix.A{j}=sqrt(mix.lambda{j}-1)'*sqrt(mix.psi(j, :)).*mix.U{j}';
            %  CM-step 3: Update \Psi
            mix.psi(j, :)=uppsi(mix, Stil, j);         
        else
            mix.lambda{j}=[];
            mix.U{j}=[];
            mix.A{j}=[];
            S1=(S(:,:,j)+S(:,:,j)')/2;
            mix.psi(j, :)=max(diag(S1)',mix.eta);
        end
    end
    % update number of parameters after reselect q
    mix.nwts = 1 + (mix.nin-1)*(mix.effdim+2) - mix.effdim.*(mix.effdim-1)/2 -1;
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
                options(8) = e; 
                disp(n);
                return;
            else  eold = e; 
            end
        end
    end
end
disp('End')
options(8) = errlog(n);
fprintf(1, 'Cycle %4d  logL %11.6f\n', n, options(8));
if display disp('Warning: Maximum number of iterations has been exceeded');end
