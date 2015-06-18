function[h]=halton(n,q);

% Based on GAUSS code from appendix of Feenberg and Skinner (ReStat, 1994)
% n= length of sequence
% q=prime number

phi=zeros(n+11,1);
i=1;

while i<n+11;
   i=i+1;
   y=1/q;
   x=1-phi(i-1);
   while x<=(y+1.e-11)
      y=y/q;
   end
   phi(i)=phi(i-1)+(q+1)*y-1;
end

h=phi(12:n+11);

clear i y x phi
