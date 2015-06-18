function r=rhaltvec(r)
%r (mxn) n draws from m bases
[m, n]=size(r);
x=rand(m, 1);
r=r+repmat(x, 1, n);
ind=r>1;
r(ind)=r(ind)-1;
