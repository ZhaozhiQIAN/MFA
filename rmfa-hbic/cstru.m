function mix1=cstru(mix1, alg)
%convert mix used in ECM to be used in EM 
if strcmp(alg, 'ECM2')
    mix1.U=[];
    mix1.effdim=[];
    mix1.lambda=[];
    Atmp=mix1.A;
    mix1.A=[];
    for j=1:mix1.ncentres
        mix1.A(:, :, j)=Atmp{j};
    end
else
    error(['Unknown algorithm']);
end