cd G:\program\my_prog\a-mfa\demo
load data.mat; % t=x;
%x=t;
options2(3)=1e-2; options2(14)=10;   options2(1)=0; options2(5)=1;%for k-means
options3(3)=1e-8; options3(14)=1000; options3(1)=0; options3(5)=0;
covar_type= 'mfa';
[ndata, xdim]=size(x); T=x'; pe=zeros(2,8);
t0=zeros(8,10); t1=t0; t2=t0;
r=1;
%cd('G:\program\my_prog\a-mfa');
%warning('off');
for ncentres=3:3  
    like1=zeros(1,10); like2=zeros(1,10); %   ncentres 
    for j=1:10       
        tic
        mix = gmmfainit1_c(x, ncentres, 2, covar_type,options2); %options=options1;        
        T0=toc; t0(ncentres,j)=T0;
        tic       
        [mix1, options, errlog] =  msqgmmfa_c(mix, T, options3,'HBIC'); % options=options3;
        T1=toc; t1(ncentres,j)=T1;
        like1(j)=options(8);    mix1tmp(j)=mix1;
%         tic
%         [mix2, options, errlog] = msqgmmfa_c(mix, T, options3,'MML');
%         T2=toc;  t2(ncentres,j)=T2;
%         like2(j)=options(8);    mix2tmp(j)=mix2;
    end
    [eopt,ind]=max(like1);    mixopt1(ncentres)=mix1tmp(ind); pe(1,ncentres)=eopt;
%    [eopt,ind]=max(like2);    mixopt2(ncentres)=mix2tmp(ind); pe(2,ncentres)=eopt;
end
% % [y,I]=max(pe,[],2); 
% % nc(3,r)=mixopt1(I(1)).ncentres; 
%nc(4,r)=mixopt2(I(2)).ncentres; 
%mixop(3,r)=mixopt1(I(1)); 
%mixop(4,r)=mixopt2(I(2));
%act=msgmmactiv(mixop(3,r), tsdata, options1); prob(3,r)=sum(log(exp(act)*(mixop(3,r).priors)'+realmin));
%act=msgmmactiv(mixop(4,r), tsdata, options1); prob(4,r)=sum(log(exp(act)*(mixop(4,r).priors)'+realmin));
T(2,r)=sum(sum(t0+t1));
%T(3,r)=sum(sum(t0+t2));


% options1(17)=1; options1(1)=0;
% options2(17)=1; options2(3)=1e-2; options2(14)=10;   options2(1)=0; options2(5)=1;%for k-means
% options3(17)=1; options3(3)=1e-6; options3(14)=1000; options3(1)=0; options3(5)=0;
%
% [ndata, xdim]=size(t);q=[[2:2:20] [21:30]];% q=[2:2:20];
% % diffs=t-repmat(mean(t,1),ndata,1); options3(21)=0.1*mean(mean(diffs.^2,2)); options2(21)=options3(21);
% for ncentres=7:25
%     ncentres
%     cd('G:\program\my_prog\imppca')
%     like=zeros(1,10);
%     for j=11:20
%         parfor k=1:10
%             options=options1;
%             % mix=gmm(10,ncentres,'ppca',5);
%             mix=msgmm(xdim,ncentres,'ppca',q(j), options);
%             options=options2;
%             % mix=gmminit(mix,t,options);
%             mix=msgmminit(mix,t,options);
%             options=options3;
%             [mix, options, errlog] = msqgmmem(mix, t, options, 'ML');
%             %[mix, options, errlog] = msqgmmaecm5(mix, t, options, 'ML');
%             % [mix, options, errlog]=gmmem(mix,t,options);
%             %disp(options(8));
%             like(k)=options(8);
%             mix1(k)=mix;
%             %         if options(8)>eopt
%             %             eopt=options(8);
%             %             mixopt1=mix;
%             %         end
%         end
%         [eopt,ind]=max(like);
%         mixopt1(ncentres,j)=mix1(ind);
%         pe(ncentres,j)=eopt-(sum(mixopt1(ncentres,j).nwts)-1)/2*log(ndata);
%
%     end
% end
% cd('G:\program\my_prog\imppca\experiment\char')
%save EM(5).mat mixopt1 pe;
%    pe
%%
%    cd('G:\program\my_prog\imppca\experiment\char')
%    load EM.mat;
%     cd('G:\program\my_prog\imppca')
%  for ncentres=7:25
%      for j=1:20
%          options=options1;    [eopt,ind]=max(like(index));
%         mixopt(ncentres)=mix1(index(ind));    post = msgmmpost(mixopt1(ncentres,j), t,options);
%         pe(ncentres,j,1)=eopt-(sum(mixopt1(ncentres,j).nwts)-1)/2*log(ndata); %BIC/MDL
%         pe(ncentres,j,2)=eopt-(sum(mixopt1(ncentres,j).nwts)-1)/2*log(ndata)+sum(log(max(post,[],2))); %ICL
%         pe(ncentres,j,3)=eopt-sum(mixopt1(ncentres,j).nwts/2.*log(max(ndata*mixopt1(ncentres,j).priors,1)))...
%         - (mixopt1(ncentres,j).ncentres-1)/2*log(ndata);% HBIC
%         pe(ncentres,j,4)=eopt-sum(mixopt1(ncentres,j).nwts/2.*log(max(ndata*mixopt1(ncentres,j).priors,1)))...
%         - mixopt1(ncentres,j).ncentres/2*log(ndata)+sum(mixopt1(ncentres,j).nwts+1)/2*(log(12)-1); % MML
%      end
%  end