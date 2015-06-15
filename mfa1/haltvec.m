function u=haltvec(m, n)
%r (mxn) n draws from m bases
u=zeros(m, n);
v=prime(1:m);
for i=1:m
    u(i, :)=halton(n, v(i));
end
