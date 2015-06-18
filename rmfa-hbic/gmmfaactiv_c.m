function [loga, Minv, psiinvA, Ax, psix2_Ax2] = gmmfaactiv_c(mix, x, x2, alg, Minv, psiinvA, Ax, psix2_Ax2)
%GMMACTIV Computes the activations of a Gaussian mixture model.
% output: Minv, psiinvA, xA to avoid re-computation in M-step.
%x, x2: dxn; Ax=qxnxM; A:qxd

ndata = size(x, 2); loga = zeros(mix.ncentres, ndata);  % Preallocate matrix
log_normal = mix.nin*log(2*pi); d2 = zeros(mix.ncentres, ndata);
logZ = zeros(mix.ncentres, 1); cal=0; if nargin<5 cal=1;end
switch alg
    case {'GEM1', 'AECM1'}
        if cal Ax=zeros(mix.subdim, ndata, mix.ncentres);%qxd
            Minv=zeros(mix.subdim, mix.subdim, mix.ncentres);%qxq
            psiinvA=zeros(mix.subdim, mix.nin, mix.ncentres);%qxdxM
            psix2_Ax2=zeros(mix.ncentres, ndata, 2);
        end
        for j = 1:mix.ncentres
            psiinv=1./mix.psi(j, :);
            if cal
                psiinvA(:, :, j)=mix.A(:, :, j).*repmat(psiinv, mix.subdim, 1);
                Minv(:, :, j)=inv(eye(mix.subdim)+mix.A(:, :, j)*psiinvA(:, :, j)');
                Ax(:, :, j)=psiinvA(:, :, j)*x;
            end
            jcentres=mix.centres(:, j);
            mupsi=jcentres'.*psiinv; %1xd
            muA=psiinvA(:, :, j)*jcentres; %qx1
            muAM=Minv(:, :, j)*muA; %qx1
            logZ(j)=log_normal+sum(log(mix.psi(j, :)))-log(det(Minv(:, :, j)))+mupsi*jcentres-muAM'*muA;
            if cal psix2_Ax2(j, :, 1)=psiinv*x2; psix2_Ax2(j, :, 2)=sum(Minv(:, :, j)*Ax(:, :, j).*Ax(:, :, j), 1);end
            d2(j, :)= psix2_Ax2(j, :, 1) + (-2*mupsi)*x - psix2_Ax2(j, :, 2)+(2*muAM)'*Ax(:, :, j);
            %d2(j, :)= psiinv*x2 + (-2*mupsi)*x-sum(Minv(:, :, j)*Ax(:, :, j).*Ax(:, :, j), 1)+(2*muAM)'*Ax(:, :, j);
        end
        loga = -0.5*(d2 + repmat(logZ, 1, ndata));
    case {'GEM2', 'ECM1', 'ECM2', 'ECM3'}
        for j = 1:mix.ncentres
            if mix.effdim(j)>0
                psiinv=1./mix.psi(j, :);
                psiinvA=mix.A{j}.*repmat(psiinv, mix.effdim(j), 1);
                Minv=inv(eye(mix.effdim(j))+mix.A{j}*psiinvA');
                Ax=psiinvA*x;
                mupsi=mix.centres(j, :).*psiinv;
                muA=psiinvA*mix.centres(j, :)';
                muAM=Minv*muA;
                logZ(j)=log_normal+sum(log(mix.psi(j, :)))-log(det(Minv))+mupsi*mix.centres(j, :)'-muAM'*muA;
                d2(j, :)= psiinv*x2 + (-2*mupsi)*x - sum(Minv*Ax.*Ax, 1)+(2*muAM)'*Ax;
            else
                psiinv=1./mix.psi(j, :);
                mupsi=mix.centres(j, :).*psiinv;
                logZ(j)=log_normal+sum(log(mix.psi(j, :)))+mupsi*mix.centres(j, :)';
                d2(j, :)= psiinv*x2 + (-2*mupsi)*x;
            end
            %            ind=find(d2(j, :)+repmat(mupsi*mix.centres(j, :)'-muAM'*muA, 1, ndata)<0);
            %             if length(ind)>0
            %                 d2(j, ind)+mupsi*mix.centres(j, :)'-muAM'*muA
            %             end
        end
        loga = -0.5*(d2 + repmat(logZ, 1, ndata));
    otherwise
        error(['Unknown algorithm ', alg]);
end


