function g=group(r)
% Group the ranking data with same pattern
% input  r: ranking data matrix N(data number)xD (dimension) 
% output g.pi: pattern index
%        g.c: pattern count
[ndata, xdim]=size(r);
D=dist2(r,r);
i=1;
ind=find(D(i,i+1:ndata)==0)+i;
g.pi=i;
g.c=length(ind)+1;
i=i+1;
while i<ndata
    if sum(i==ind)==0
        ind1=find(D(i,i+1:ndata)==0)+i;
        ind=[ind ind1];
        g.pi=[g.pi i];
        g.c=[g.c; length(ind1)+1];
    end
    i=i+1;
end
if sum(i==ind)==0
    g.pi=[g.pi i];
    g.c=[g.c; 1];
end
