function [WC_S,WC_M,WC_L,TC_S,TC_M,TC_L,TF_C,...
    ND_S,ND_M,ND_L,TN_C,CT_S,CT_M,CT_L,CT_T,CR_S,CR_M,CR_L,CR_T] = dynamic_networks(data,nsim,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates Dyanmic Networks introduced in:                 %
% Barunik J. and Ellington M. (2020): Dynamic Networks in Large Financial %
% and Economic Systems                                                    %
%                                                                         %
% and                                                                     %
%                                                                         %
% Barunik J. and Ellington M. (2020): Dynamic Network Risk                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All code developed by Michael Ellington, Jozef Barunik and Lubos Hanus  %
% This was written August 26, 2020 By Michael Ellington                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data: data matrix                                                       %
% nsim: Number of simulations                                             %
% L:    lag length                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WC_S: Within Short-term connectedness
% WC_M: Within Medium-term connectedness
% WC_L: Within Long-term connectedness
% TC_S: Short-term connectedness
% TC_M: Medium-term connectedness
% TC_L: Long-term connectedness
% TF_C: Total Frequency connectedness
% ND_S: Short-term NET directional connectedness
% ND_M: Medium-term NET directional connectedness
% ND_L: Long-term NET directional connectedness
% TN_C: Total NET directional connectedness
% CT_S: Short-term TO connectedness
% CT_M: Medium-term TO connectedness
% CT_L: Long-term TO connectedness
% CT_T: Total TO connectedness
% CR_S: Short-term FROM connectedness
% CR_M: Medium-term FROM connectedness
% CR_L: Long-term FROM connectedness
% CR_T: Total FROM connectedness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE THAT IN THE FUNCTION get_dynnet() WE DEFINE FREQUENCIES BASED ON OUR
% APPLICATION OF DAILY DATA. YOU CAN MODIFY ACCORDINGLY TO EXAMINE ANY
% HORIZON OF INTEREST.

% In our case: short-term is defined as 1-day to 1-week
%              medium-term is 1-week to 1-month
%              long-term is horizons > 1-month

rng('default'); rng(1);
[T,N]=size(data);
shrinkage=0.05; % Overall shrinkage of Minnesota Prior.
% Get data in VAR setup
dat2=data;
for i=1:L
    temp=lag0(dat2,i);
    X(:,1+N*(i-1):i*N)=temp(1+L:T,:);
end
y=dat2(1+L:T,:); T=T-L; X=[ones(T,1),X]; % y[T x N], X[T x N*L+1]
clear temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SI,PI,a,RI]=Minn_NWprior(dat2,T,N,L,shrinkage);
% SI prior mean matrix. PI prior variance of SI (diagonal matrix (K x K)),
% and RI is prior for covariance matrix of VAR model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate weights
%weights=normker(T,sqrt(T)); % This is PETROVA'S DEFAULT
weights=normker(T,8); % HERE 8 Denotes the width of the kernel change if you wish
% Follow Petrova (2019) work with precision
priorprec0=PI^(-1);
clear PI data dat2 ind
priorprec0=sparse(priorprec0); % create sparse matrix allowing more efficient allocation of memory.
RI=sparse(RI); % create sparse matrix allowing more efficient allocation of memory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posteriors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stab_ind=zeros(nsim,T,'single');
max_eig=zeros(nsim,T,'single');
% Generate storage matrices
HO=10+1; % DEFAULT FOR PACKAGE ONLY. WHEN USING IN RESEARCH CHANGE TO HO=100+1
         % HOWEVER PLEASE NOTE THAT THIS WILL TAKE MORE TIME
diagonal=0; % set to 0 for full covariance matrix, 1 for diagonal.
            % In our example we set diagonal to zero. However as dimension
            % of VAR gets large Total Connectedness tends to 1
TC_S=single(zeros(T,nsim));
TC_M=single(zeros(T,nsim));
TC_L=single(zeros(T,nsim));
WC_S=TC_S; WC_M=TC_M; WC_L=TC_L; % Within connectedness
TF_C=single(zeros(T,nsim));
ND_L=single(zeros(T,N,nsim));
ND_M=ND_L; ND_S=ND_L; TN_C=ND_L; % Net-directional connectedness over frequency bands.
CT_L=ND_L; CT_M=ND_L; CT_S=ND_L; CT_T=ND_L;
CR_L=ND_L; CR_M=ND_L; CR_S=ND_L; CR_T=ND_L;

tic;

parfor kk=1:T
%   kk/T % This reports the proportion of computations carried out.
   % when using parallel computing numbers will not be in order.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Storage for connectedness measures
   wcs=zeros(1,nsim);
   wcm=wcs; wcl=wcs; 
   tcs=wcs; tcm=wcs; tcl=wcs; 
   tfc=zeros(1,nsim); 
   ndcl=zeros(N,nsim); ndcm=ndcl; ndcs=ndcl; tndc=ndcl;
   ctcl=zeros(N,nsim); ctcm=ndcl; ctcs=ndcl; ctct=ndcl;
   crcl=zeros(N,nsim); crcm=ndcl; crcs=ndcl; crct=ndcl;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Estimation and connectedness computing
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   w=weights(kk,:);
   bayesprec=(priorprec0+X'*diag(w)*X);
   bayessv=bayesprec^(-1);
   BB=bayessv*((X'*diag(w))*y+priorprec0*SI);
   bayesb=BB(:);
  bayesalpha=a+sum(w);
  g1=SI'*priorprec0*SI;
  g2=y'*diag(w)*y;
  g3=BB'*bayesprec*BB;
  bayesgamma=RI+g1+g2-g3;
  bayesgamma=0.5*bayesgamma+0.5*bayesgamma'; % it is symmetric but just in case
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Draws
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ii=1:nsim
     mm=0;
     while mm<1
         SIGMA=iwishrnd(bayesgamma,bayesalpha); % Draw from IW distribution
         nu=randn(N*L+1,N);
         Fi1=(BB+chol(bayessv)'*nu*(chol(SIGMA)))';         
         max_eig(ii,kk)=max(abs(eig([Fi1(:,2:end); eye(N), zeros(N,N)])));
         if max_eig(ii,kk)<.999 % check stability of draw
             stab_ind(ii,kk)=1;
             mm=1;
         end
     end
     [~,wold]=get_GIRF(Fi1,SIGMA,1,L,HO-1);
     [ttfc,tc1,tc2,tc3,wc1,wc2,wc3,ct1,ct2,ct3,ctt,cr1,cr2,cr3,crt,... 
         ndc1,ndc2,ndc3,tnd1]=get_dynnet(wold,T,SIGMA,diagonal); % efficient code to get frequency connectedness
     wcs(:,ii)=wc3; wcm(:,ii)=wc2; wcl(:,ii)=wc1;
     tcs(:,ii)=tc3; tcm(:,ii)=tc2; tcl(:,ii)=tc1;
     tfc(:,ii)=ttfc;
     ndcl(:,ii)=ndc1; ndcm(:,ii)=ndc2; ndcs(:,ii)=ndc3; tndc(:,ii)=tnd1;
     ctcl(:,ii)=ct1; ctcm(:,ii)=ct2; ctcs(:,ii)=ct3; ctct(:,ii)=ctt;
     crcl(:,ii)=cr1; crcm(:,ii)=cr2; crcs(:,ii)=cr3; crct(:,ii)=crt;
  end
 WC_S(kk,:)=wcs; WC_M(kk,:)=wcm; WC_L(kk,:)=wcl;
 TC_S(kk,:)=tcs; TC_M(kk,:)=tcm; TC_L(kk,:)=tcl;
 TF_C(kk,:)=tfc; 
 ND_S(kk,:,:)=ndcs; ND_M(kk,:,:)=ndcm; ND_L(kk,:,:)=ndcl;
 TN_C(kk,:,:)=tndc;
 CT_S(kk,:,:)=ctcs; CT_M(kk,:,:)=ctcm; CT_L(kk,:,:)=ctcl;
 CT_T(kk,:,:)=ctct;
 CR_S(kk,:,:)=crcs; CR_M(kk,:,:)=crcm; CR_L(kk,:,:)=crcl;
 CR_T(kk,:,:)=crct;
end
toc
end