function [y,p,tt, wgt] = ghkr(mu,varcov,u,I)
%GHK routine

%	INPUTS:	mu (dx1)= mean column vector of mvn variable vector
%		 sigma (dxd)= var/cov matrix
%		     u (d-1x r)= uniform rand number 
%			I  = initialization of parameters
%				
%
%	OUTPUT:	p = ranking prob.
%        y.y1 = E(x); 
%        y.y2 = E(xx'); 
%
%	Originally written by:
%	Vassilis A. Hajivassiliou
%	Source: 'Some Practical Issues in Maximum Simulated Likelihood',
%	Working Paper, London School of Economics, October 1997, p.17.
%
%	This version written by:
%	Jinahua Zhao
%	Department of statistics & actuarial science
%	The University of Hong Kong
%	jhzhao.ynu@gmail.com

m=length(mu);r=size(u, 2);%r = number of draws
c=chol(varcov, 'lower');  %get lower triangular Cholesky for varcov
mint=7.6506;minc=1e-14;upt=mint*I.e; tt=I.t;

j=1;
bi=-mu(j)/c(j,j);
tb=1;
if bi<-9  %only compute \phi(x) in -9<x<9 to reduce evaluation of \phi
    tb=0;
    tt(1, :)=-upt;
else if bi<9
        tb=phi(bi);
    end
    ac=u(j,:)*tb;t=upt; 
    ind3=ac<minc;    t(ind3)=-mint;
    ind4=(ac>=minc & ac<1-minc);% only compute \phi^{-1}(y) in -minc<y<1-minc to reduce evaluation of \phi^{-1}
    t(ind4)=phinv(ac(ind4));
    tt(1, :)=t;
end
wgt=tb;
while j<m-1
    j=j+1;tb=I.e;t=upt;
    bi=(-mu(j)-c(j,1:j-1)*tt(1:j-1, :))/c(j,j);
    ind1=bi<-9;    tb(ind1)=0;
    ind2=(bi>=-9 & bi<9);    tb(ind2)=phi(bi(ind2)); 
    ac=u(j,:).*tb;     ind3=ac<minc;      t(ind3)=-mint;
    ind4=(ac>=minc & ac<1-minc);
    t(ind4)=phinv(ac(ind4));     
    tt(j, :)=t;
    wgt=wgt.*tb;
    %    end
end

%sum 
p=sum(wgt);
ct=c(1:m-1, 1:m-1);
%compuate E(x_{d-1}) and E(x_{d-1}x_{d-1}')
y.y1=ct*(tt*wgt')/p;
tt2=tt.*repmat(sqrt(wgt),m-1,1);
y.y2=ct*(tt2*tt2'/p)*ct';
%compuate E(x) and E(xx')
sigma11=varcov(1:m-1, 1:m-1);
sigma12_1=varcov(m, 1:m-1)*inv(sigma11);
y.y1=[y.y1;sigma12_1*y.y1];
y12=y.y2*sigma12_1';
y22=varcov(m,m)-sigma12_1*(sigma11-y.y2)*sigma12_1';
y.y2(1:m-1, m)=y12;
y.y2(m, 1:m-1)=y12';
y.y2(m,m)=y22;
p=p/r;
return

function p = phi(z)
%
%  Standard statistical normal distribution
%
p = erfc( -z/sqrt(2) )/2;
return
%
% end phi
%
function z = phinv(w)
%
%  Standard statistical inverse normal distribution
%
z = -sqrt(2)*erfcinv( 2*w );
return
%
% end phinv