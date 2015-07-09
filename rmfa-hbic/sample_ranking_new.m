%% generate sample
dim=15;
ndata=2000;
% change mfa.subdim into a vector
dppca_dim=[2,2];
mfa.ncentres=2;      
% change mfa.A from 3d array into list of 2d-arrays
mfa.psi=zeros(mfa.ncentres, dim);    
mfa.centres=0.2*[1:dim;dim:-1:1]';
mfa.A{1} = [repmat([2.8,0.8;2.8,-0.8],7,1);[2.8,0.8]];
mfa.A{2} = mfa.A{1};

for j=1: mfa.ncentres
    mfa.psi(j, :)=j*rand(1, dim);
end
% mfasamp: sample data from mfa model
mfa.subdim=dppca_dim; 
x=mfasamp(mfa, ndata);
%% generate rankings
[y, ro]=sort(x, 1, 'ascend');
g=group(ro');
% r stores rankings: r.p is ranking format (column wise), r.c is count vector 
r.p=ro(:, g.pi);
r.c=g.c;