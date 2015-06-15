function Psi=uppsi(mix, B, k)
%UPPSI update \Psi using the improved algorithm described in [1].
%
%Input: mix, \btS
%Output: mix with new mix.psi.


%   reference: [1] ML Estimation for Factor Analysis: EM or non-EM? by Jianhua Zhao,
%   Philip L. H. Yu and Qibao Jiang. Department of Statistics and Actuarial
%   Science, The University of Hong Kong

%	Copyright (c) Jianhua Zhao (2007)

xdim=mix.nin;
etat=mix.eta./mix.psi(k, :)-1;
C=eye(xdim); lambda_inv=1./mix.lambda{k}; U=mix.U{k}.*repmat(lambda_inv-1, xdim, 1);
V=mix.U{k}';  lambda_diff=lambda_inv-mix.lambda{k}; T=zeros(xdim, 1); %Brow=zeros(1, xdim);
for j=1: floor(xdim/2)
    Crow=C(j, 1:j);
    Vcol=V(:, 1:j)* Crow';
    T(j:xdim)=U(j:xdim, :)*Vcol;
    a=T(j)+1;
    b=lambda_diff.*Vcol'*Vcol+ Crow*B(1:j, 1:j)*Crow';
    %Brow(1:j)=Crow(1:j);b=lambda_diff.*Vcol'*Vcol+ Brow*B*Brow'; %much less
    %efficient than line 22.
    omega=max((b-a)/a^2, etat(j));
    if omega>etat(j) || mix.psi(k, j)~=mix.eta % if old \psi=eta and \omega_i=etat_i, no need to update \psi and C
        mix.psi(k, j)=(omega+1).*mix.psi(k, j);
        r=omega/(1+omega*a);
        C(j+1:xdim, 1:j)= C(j+1:xdim, 1:j)-r*T(j+1:xdim)* Crow;
    end
end
C=C'; Bcol=zeros(xdim, 1);
for j=floor(xdim/2)+1:xdim
    Ccol=C(1:j, j);
    Vcol=V(:, 1:j)* Ccol;
    T(j:xdim)=U(j:xdim, :)*Vcol;
    a=T(j)+1;    Bcol(1:j)=Ccol(1:j);b=lambda_diff.*Vcol'*Vcol+ Bcol'*B*Bcol;  %
    %b=lambda_diff.*Vcol'*Vcol+ Ccol'*B(1:j, 1:j)*Ccol; %less
    %efficient than line 37.
    omega=max((b-a)/a^2, etat(j));
    if omega>etat(j) || mix.psi(k, j)~=mix.eta % if old \psi=eta and \omega_i=etat_i, no need to update \psi and C
        mix.psi(k, j)=(omega+1).*mix.psi(k, j);
        if j<xdim
            r=omega/(1+omega*a);
            C(1:j, j+1:xdim) = C(1:j, j+1:xdim)-r*Ccol*T(j+1:xdim)';
        end
    end
end
Psi=mix.psi(k, :);
