%% DynamicNets 
% Master file providing example of how to estimate dynamic networks
% introduced in:

% Barunik J. and Ellington M. (2020): Dynamic Networks in Large Financial
% and Economic Systems

% and

% Barunik J. and Ellington M. (2020): Dynamic Network Risk

% We use the TVP VAR proposed by Katerina Petrova in:

% Petrova, K. (2019) A quasi-Bayesian local likelihood approach to 
% time varying parameter VAR models. Journal of Econometrics 212, 286--306

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All code developed by Michael Ellington, Jozef Barunik and Lubos Hanus  %
% This was written August 26, 2020 By Michael Ellington                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; 
addpath('functions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data and Preliminaries                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('exampledata.txt')
data=exampledata(:,1:4); % N=4 variable VAR model
clear exampledata
data=sqrt(data); % data matrix
nsim=100;       % number of simulations
L=2;             % lag length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Dynamic Network                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[WC_S,WC_M,WC_L,TC_S,TC_M,TC_L,TF_C,...
 ND_S,ND_M,ND_L,TN_C,CT_S,CT_M,CT_L,...
 CT_T,CR_S,CR_M,CR_L,CR_T]=dynamic_networks(data,nsim,L);
% Time to estimate on 64GB Desktop PC with i7 3.70GHz 6-core processor:
% 226.22 seconds.

% NOTE HERE THAT WE HAVE SOME DEFAULTS IN THE SCRIPT dynamic_networks
% Note these outputs are the full posterior distribution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now define 95% confidence intervals of the posterior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qq=[0.025 0.5 0.975];
WC_S=quantile(WC_S,qq,2);
WC_M=quantile(WC_M,qq,2);
WC_L=quantile(WC_L,qq,2); % WITHIN CONNECTEDNESS
TC_S=quantile(TC_S,qq,2);
TC_M=quantile(TC_M,qq,2);
TC_L=quantile(TC_L,qq,2); 
TF_C=quantile(TF_C,qq,2); % FREQUENCY CONNECT
ND_S=quantile(ND_S,qq,2);
ND_M=quantile(ND_M,qq,2);
ND_L=quantile(ND_L,qq,2); 
TN_C=quantile(TN_C,qq,2); % NET-DIRECTIONAL FREQUENCY CONNECT
CT_S=quantile(CT_S,qq,3);
CT_M=quantile(CT_M,qq,3);
CT_L=quantile(CT_L,qq,3); 
CT_T=quantile(CT_T,qq,3); % TO FREQUENCY CONNECT
CR_S=quantile(CR_S,qq,3);
CR_M=quantile(CR_M,qq,3);
CR_L=quantile(CR_L,qq,3); 
CR_T=quantile(CR_T,qq,3); % FROM FREQUENCY CONNECT

% _S denotes short-term, _M denotes medium-term, _L denotes long-term
% _T/_C denotes total (i.e. sums over frequency bands)

% In this example:
% Short-term is 1 day to 5-days
% Medium-term is 5-days to 20-days
% Long-term is > 20-days

% Depending on application you will need to change in get_dynnet() in the 
% functions folder 

% If you want full posterior distribution delete/comment out lines 46-65
% However if you do so the plots below will not work!
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot some results                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=(1:1:length(data)-L)';
figure(1)
plotx4(t,[TF_C(:,2), TF_C(:,1), TF_C(:,3)])
axis tight
ylim([0 100])
title('Total Frequency Connectedness, C_{\infty}')

figure(2)
plot(t,TF_C(:,2),'k','LineWidth',1.3)
hold on,
plot(t, TC_S(:,2),'color',[0 0.5 0],'LineWidth',1.3)
plot(t, TC_M(:,2),'b','LineWidth',1.3)
plot(t, TC_L(:,2),'r','LineWidth',1.3)
plotx4(t,[TF_C(:,2), TF_C(:,1), TF_C(:,3)])
plotx3(t,[TC_S(:,2), TC_S(:,1), TC_S(:,3)])
plotx6(t,[TC_M(:,2), TC_M(:,1), TC_M(:,3)])
plotx5(t,[TC_L(:,2), TC_L(:,1), TC_L(:,3)])
axis tight
ylim([0 100])
title('TVP Network Connectedness')
legend('Total','Short','Medium','Long','Location','SouthOutside')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save outputs                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFILE='Results_DynamicNets';
varname(1,:)='WC_S'; varname(2,:)='WC_M'; varname(3,:)='WC_L';
varname(4,:)='TC_S'; varname(5,:)='TC_M'; varname(6,:)='TC_L';
varname(7,:)='TF_C'; varname(8,:)='ND_S'; varname(9,:)='ND_M'; 
varname(10,:)='ND_L'; varname(11,:)='TN_C'; varname(12,:)='CT_S';
varname(13,:)='CT_M'; varname(14,:)='CT_L'; varname(15,:)='CT_T';
varname(16,:)='CR_S'; varname(17,:)='CR_M'; varname(18,:)='CR_L'; 
varname(19,:)='CR_T';

save(DFILE,varname(1,:),varname(2,:),varname(3,:),varname(4,:),varname(5,:),...
    varname(6,:),varname(7,:),varname(8,:),varname(9,:),varname(10,:),varname(11,:),...
    varname(12,:),varname(13,:),varname(14,:),varname(15,:),varname(16,:),...
    varname(17,:),varname(18,:),varname(19,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
