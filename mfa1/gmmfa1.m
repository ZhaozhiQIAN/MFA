function [mix, options, errlog, ndraw] = gmmfa1(mix, r, options, alg)
%GMMFA	EM algorithms for Gaussian mixture model.
%increase no. of draws if L^(t+1)<L^t
%x, x2: dxn; Ax=qxnxM; A:qxd; psiinvA:qxd; loga: MxN; post: MxN;
%options(15): number of draw;

% Sort out the options
ndrawi=options(21);ndraw_max=options(22);
ndata=sum(r.c);[xdim, npat]= size(r.p);

if options(14) niters = options(14);else niters = 100;end
display = options(1);store = 0;
if (nargout > 2) store = 1;	% Store the error values to return them
    errlog = zeros(1, niters);
    ndraw = errlog;
end
ndraw(1:2)=options(20);
test = 0;
if options(3)
    test = 1;
end

if (options(14))
    niters = options(14);
else
    niters = 100;
end
C0=zeros(xdim);I.C=C0;
for j=1:xdim-1
    I.C(j, j: j+1)=[1 -1];
end
I.C(xdim, xdim)=1;I.Ci=inv(I.C);I.C0=zeros(xdim, xdim, mix.ncentres);
I.p=zeros(mix.ncentres,npat);I.x1_t=zeros(xdim,mix.ncentres);
I.e=ones(1, ndraw_max);I.t=zeros(xdim-1, ndraw_max);

u=haltvec(xdim-1, ndraw_max);

% Main loop of algorithm
for n = 1:niters
    options(16)=n;
    % calling function estep1
    [mix, S, e]=estep1(r, mix, I, u(:, 1:ndraw(n)), alg);
    if (display || store || test)
        if store  errlog(n) = e;
            if mod(n,10)==0
                fprintf(1, 'Cycle %4d  logL %11.6f\n', n, e);
            end
        end
        if (n > 1 & e<eold)
            ndraw(n+1)=min(ndraw(n)+ndrawi, ndraw_max);
        else
            ndraw(n+1)=ndraw(n);
        end
        if test
            if (n > 1 & abs(e - eold) < options(3))
                options(8) = e; return;
            else
                eold=e;
            end
        end
    end
    switch alg
        case {'ECM1','ECM2','ECM3'}
            for j=1: mix.ncentres
                %CM-step 2; Update A
                psisinv=mix.psi(j, :).^-0.5;
                St=psisinv'*psisinv.*S(:, :, j);
                [templambda, tempU] = eigdec(St, mix.subdim);
                mix.lambda{j}=templambda(templambda>1)';
                mix.effdim(j)=length(mix.lambda{j});
                mix.U{j}=tempU(:, 1:mix.effdim(j));
                mix.A{j}=sqrt(mix.psi(j, :))'*sqrt(mix.lambda{j}-1).*mix.U{j};
                %sfa.A=repmat(sqrt(sfa.psi)', 1, sfa.effdim).*(sfa.U.*repmat(sqrt(sfa.lambda-1), xdim, 1));

                switch alg
                    case 'ECM1'
                        mix.psi(j, :)=duppsi(mix, St, j, alg);    %CM-step 3: Update \Psi
                    case 'ECM2'
                        mix.psi(j, :)=uppsi(mix, St, j);         %CM-step 3: Update \Psi
                    case 'ECM3'
                        mix.psi(j, :)=luppsi(mix, diag(St), j, alg);         %CM-step 3: Update \Psi
                end
            end
        % relavent case
        case 'EM'
            for j=1: mix.ncentres
                psiinv=mix.psi(j, :).^-1;
                psiinvA=repmat(psiinv', 1, mix.subdim).*mix.A(:, :, j);
                Minv=inv(eye(mix.subdim)+mix.A(:, :, j)'*psiinvA);
                SpsiinvA=S(:, :, j)*psiinvA;
                SpsiinvAMinv=SpsiinvA*Minv;
                mix.A(:, :, j)=SpsiinvA*inv(eye(mix.subdim)+SpsiinvAMinv'*psiinvA);
                mix.psi(j, :)=max(diag(S(:, :, j))-sum(SpsiinvAMinv.*mix.A(:, :, j), 2), mix.eta)';
            end
    end
end
options(8) = errlog(n);
%fprintf(1, 'Cycle %4d  logL %11.6f\n', n, options(8));

