function [TFC,TC1,TC2,TC3,WC1,WC2,WC3,DCTL,DCTM,DCTS,DCTT,DCRL,DCRM,DCRS,DCRT,...
    NDCL,NDCM,NDCS,TNDC]=get_dynnet(wo,TT,sig,diagonal)
%**************************************************************************
% Michael Ellington 26/08/2020
% Dynamic Network Estimation as in:
% Barunik J. and Ellington M. (2020): Dynamic Networks in Large Financial
% and Economic Systems
%
% and
%
% Barunik J. and Ellington M. (2020): Dynamic Network Risk
%
%
% Inputs: wo is wold decomposition of MA coefficients [N x N x HO]
%
%         TT is the number of time series observations used in estimating
%         TVP VAR model.
%
%         sig is a draw of the estimated posterior distribution of the
%         covariance matrix of residuals stemming from the TVP VAR at time
%         t. [N x N]
%
%         diagonal is binary: 0 means use full covariance matrix. 1 means
%         take diagonal.
%
% Outputs: TFC is total frequency connectedness
%          TC1 is long-term connectedness
%          TC2 is med-term connectedness
%          TC3 is short-term connectedness
%          WC1 is within long-term connectedness
%          WC2 is within med-term connectedness
%          WC3 is within short-term connectedness
%          DCTL is long-term to connectedness
%          DCTM is med-term to connectedness
%          DCTS is short-term to connectedness
%          DCTT is total net to connectedness
%          DCRL is long-term from connectedness
%          DCRM is med-term from connectedness
%          DCRS is short-term from connectedness
%          DCRT is total net from connectedness
%          NDCL is long-term net directional connectedness
%          NDCM is med-term net directional connectedness
%          NDCS is short-term net directional connectedness
%          TNDC is total net directional freq connectedness
%**************************************************************************
% Define frequency window;
Tw=floor(TT/10); % approximates the total frequency window.
%Tw=10;
omeg=linspace(0,pi,Tw)'; % create equally spaced line from 0 to pi in Tw intervals.

% Define bands
omeg2=pi./omeg; 
bandmat=[omeg,omeg2];
d1=bandmat(bandmat(:,2)>20); % long term equals (20,260+] days
d2=bandmat(bandmat(:,2)<=20 & bandmat(:,2)>5); % medium term equals (5,20] days
d3=bandmat(bandmat(:,2)<=5); % short term equals [1,5] days
% d1 and d2 are not strictly open, since double counting when summing over
% frequencies will violate condition in line 165.

d1=length(d1); d2=length(d2); d3=length(d3);

% If omeg is frequency window, we decompose into three components:
% To get long-term we sum from 1:length(d1)
% To get med-term we sum from length(d1)+1:length(omeg)-length(d2)
% To get short-term we sum from
% length(d1)+length(d2)+1:end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here instead of frequency band, we are looking at (w={0,...,pi}) Which,
% due to symmetry of integral (-pi, pi) is the whole window...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HO=size(wo,3); % 3rd dimension of wold coefficients (i.e. HOR)
N=size(wo,1);  % Number of variables (and shocks) in model
i=sqrt(-1);

% First get Omega = \sum_{w} \Psi(w)Sigma\Psi'(w) 
% and also define Fourier transforms of coefficient matrices (wold
% decomposition).
Omeg=zeros(N,N,length(omeg));
GI=zeros(N,N,length(omeg));
for w = 1:length(omeg)
    for nn=1:HO
       GI(:,:,w)=GI(:,:,w)+wo(:,:,nn)*exp(-i*nn*omeg(w));  
    end
    if diagonal==0
    Omeg(:,:,w)=GI(:,:,w)*sig*GI(:,:,w)'; 
    elseif diagonal==1
    Omeg(:,:,w)=GI(:,:,w)*diag(diag(sig))*GI(:,:,w)'; 
    end
end
Omeg=sum(real(Omeg),3); % This is overall omega (spec density over whole frequency window)
% The diagonal elements of this matrix are (Omeg)j,j that we use in line
% 96.


% First define numerator over all frequency windows
FC=zeros(N,N,length(omeg)); %FC1=FC;
for w=1:length(omeg)
   if diagonal==0
   GI1=GI(:,:,w)*sig; % product of GI(w) and Sig at given frequency
   elseif diagonal==1
   GI1=GI(:,:,w)*diag(diag(sig)); 
   end
   PS=zeros(N,N); PP=zeros(N,N);
   for k=1:N % following Defintion 2.3 in Barunik and Krehlik (2018)
       for j=1:N % jth row and kth column
        PS(j,k)=(abs(GI1(j,k)))^2;
        PP(j,k)=PS(j,k)/(Omeg(j,j)*sig(k,k)); % Use Omeg as defined earlier and take ratio.
       end
   end  
FC(:,:,w)=PP;
end

PP1=sum(FC,3);
for w=1:length(omeg)
for j=1:N
   FC(j,:,w)=FC(j,:,w)./sum(PP1(j,:)); % normalise theta by row sum of theta_inf 
end
end

thetainf=sum(FC,3);

% Now get long-term
temp1=sum(FC(:,:,1:d1),3); % theta_{d1} summed over (0 to 0.1576)

% get net directional connectedness matrix on Frequency L
for i=1:N
   tr(i,:)=sum(temp1(i,:))-temp1(i,i);
   tt(i,:)=sum(temp1(:,i))-temp1(i,i);
   DCRL(i,:)=sum(temp1(i,:));
   DCTL(i,:)=sum(temp1(:,i));
%   DCRL(i,:)=sum(temp1(i,:))/temp1(i,i);
%   DCTL(i,:)=sum(temp1(:,i))/temp1(i,i); % normalises directions to own shocks
end
NDCL=(tt-tr);

WC1=100*(1-trace(temp1)/sum(sum(temp1)));
TC1=WC1*(sum(sum(temp1))/sum(sum(thetainf)));

% Now get med-term % theta_{d2} summed over (0.1576 to 0.6283)
temp1=sum(FC(:,:,d1+1:length(omeg)-d3),3);

% get net directional connectedness matrix on Frequency M
for i=1:N
   tr(i,:)=sum(temp1(i,:))-temp1(i,i);
   tt(i,:)=sum(temp1(:,i))-temp1(i,i);
   DCRM(i,:)=sum(temp1(i,:));
   DCTM(i,:)=sum(temp1(:,i));
%   DCRM(i,:)=sum(temp1(i,:))/temp1(i,i);
%   DCTM(i,:)=sum(temp1(:,i))/temp1(i,i);
end
NDCM=(tt-tr);

WC2=100*(1-trace(temp1)/sum(sum(temp1)));
TC2=WC2*(sum(sum(temp1))/sum(sum(thetainf)));

% Now get short-term % theta_{d3} summed over (0.6283 to \pi)
temp1=sum(FC(:,:,d1+d2+1:end),3);

% get net directional connectedness matrix on Frequency S
for i=1:N
   tr(i,:)=sum(temp1(i,:))-temp1(i,i);
   tt(i,:)=sum(temp1(:,i))-temp1(i,i);
   DCRS(i,:)=sum(temp1(i,:));
   DCTS(i,:)=sum(temp1(:,i));
%   DCRS(i,:)=sum(temp1(i,:))/temp1(i,i);
%   DCTS(i,:)=sum(temp1(:,i))/temp1(i,i);
end
NDCS=(tt-tr);

WC3=100*(1-trace(temp1)/sum(sum(temp1)));
TC3=WC3*(sum(sum(temp1))/sum(sum(thetainf)));

% Total Frequency Connect
tfc=TC1+TC2+TC3;

% below is a checker that proposition holds.
temp1=sum(FC,3);
% get net directional connectedness matrix on Frequency Total
for i=1:N
   tr(i,:)=sum(temp1(i,:))-temp1(i,i);
   tt(i,:)=sum(temp1(:,i))-temp1(i,i);
   DCRT(i,:)=sum(temp1(i,:));
   DCTT(i,:)=sum(temp1(:,i));
%   DCRT(i,:)=sum(temp1(i,:))/temp1(i,i);
%   DCTT(i,:)=sum(temp1(:,i))/temp1(i,i);
end
TNDC=(tt-tr);

% Total Frequency Connect: should be same as line 142.   
 tfc1=100*(sum(sum(temp1))/sum(sum(thetainf))-trace(temp1)/sum(sum(thetainf)));

if abs(tfc-tfc1)<eps+1.0e-10 % If these are equal to approximately 10d.p then allow function to go through.
   TFC=tfc;
else
    error('DIFFERENCE EXCEEDS ALLOWANCE OF eps+1.0e-10')
end

end