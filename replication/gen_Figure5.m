% -------------------------------------------------------------------------
% Plot Figure 5 in the Paper (Structural Shocks e)
% -------------------------------------------------------------------------

restoredefaultpath
clear; clc; close all;
addpath './supporting_code'

%% Load Data
IRF_lim = 60; % IRF length
lag = 6;      % VAR Lag

% Load Data
U_data = importdata('Replication_Data.xlsx');
U_data = U_data.data.Data;
dates = 1960 + 7/12 : 1/12 : 2015+4/12;
Uf = U_data(:,1);
Um = U_data(:,2);
ip = U_data(:,4);
Vs = U_data(:,6);
Gold_o = U_data(:,7);
Gold_o = [NaN;diff(log(Gold_o))];

%% Trim Sample (If needed)
ind = isnan(Um) ==0;
Um = Um(ind);
Uf = Uf(ind);
ip = ip(ind);
Vs = Vs(ind);
Gold_o = Gold_o(ind);
dates = dates(ind);
dates_lag = dates(lag+1:end);

%% Stack variables
X = [Um,ip,Uf];
X_var = X(lag+1:end,:);
X_exp = [];
for jjj = 1:lag
    X_exp= [X_exp, X(lag+1-jjj:end-jjj,:)];
end

%% Reduced Form VAR
T = size(X_var,1);
data = [Vs(lag+1:end), Gold_o(lag+1:end)];
reg31 = regstats(X_var(:,1), X_exp);
reg32 = regstats(X_var(:,2), X_exp);
reg33 = regstats(X_var(:,3), X_exp);
eta_m_OLS = reg31.r;
eta_y_OLS = reg32.r;
eta_f_OLS = reg33.r;
eta = [eta_m_OLS,eta_y_OLS,eta_f_OLS];
eta_vec_OLS = eta;
vec_A = [reg31.beta(2:end);reg32.beta(2:end);reg33.beta(2:end)];
const_var = [reg31.beta(1), reg32.beta(1), reg33.beta(1)];

%% Generate Random Rotation
load QQ_full % Load pre-saved seeds used in the paper

% Comment out to generate new random rotations
% QQ =[];
% NQ= 1.5*10^(6); %Number of Random Rotation
% parfor iii = 1:NQ
%     %generate a random rotation
%     v = randn(3,3); %First Generate a totally random 3x3 matrix v
%     [q, r]=qr(v,0); %Obtain the orthonormal matrix of v, call it q
%     QQ(:,:,iii)=q*diag(sign(diag(r)));    % this makes the diagonal of r positive
% end
%% Compute e
Ebar = [4.1634,4.0475,4.5672,4.7314];
% Generate B
B_out= gen_B(vec_A, QQ, X_var, X_exp,const_var,data,dates_lag,lag,Ebar);
% Gen ehat
[~,~,~,B_allsol,~,ehat_maxG]=gen_lo_hi(vec_A, B_out,X_var, X_exp,const_var,data,lag);

%% Left Panel

fig = plot_ts_eshock(ehat_maxG',dates(1+lag:end),[1,1,1]');
set(gcf,'paperpositionmode','manual','paperunits','inches');
dim = [10,8];
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure5_left');

%% Right Panel

B_list=  B_allsol;
dates = (1960 + 7/12 : 1/12 : 2015+4/12)'; 
dates_lag = dates(lag+1:end);
 
 
Nq = size(B_list,3);
e_list = [];
for ii = 1:Nq
e_list(:,:,ii) = (inv(B_list(:,:,ii))*eta')';
end
dates = dates_lag;
fig = figure;
subplot(3,1,1)
 
for ii = 1:Nq
plot(dates, e_list(:,1,ii), 'b', 'Linewidth',3)
hold on
end
plot(dates, 3*ones(size(dates,1),1), 'k--', 'Linewidth',3)
title('Positive $e_{\mathrm{M}}$ exceeding 2 standard deviations','interpreter','latex','Fontsize',15)
max_e1 = max(max(e_list(:,1,:)))+0.01;
ylim([2,5.55])
xlim([dates(1),dates(end)]);
rshade(dates);
subplot(3,1,2)
for ii = 1:Nq
plot(dates, -e_list(:,2,ii), 'b', 'Linewidth',3)
hold on
end
plot(dates, 3*ones(size(dates,1),1), 'k--', 'Linewidth',3)
title('Negative $e_{\mathrm{ip}}$ exceeding 2 standard deviations','interpreter','latex','Fontsize',15)
max_e2= min(min(e_list(:,2,:)))-0.01;
ylim([2,5.5])
xlim([dates(1),dates(end)]);
rshade(dates);
subplot(3,1,3)
for ii = 1:Nq
plot(dates, e_list(:,3,ii), 'b', 'Linewidth',3)
hold on
end
plot(dates, 3*ones(size(dates,1),1), 'k--', 'Linewidth',3)
title('Positive $e_{\mathrm{F}}$ exceeding 2 standard deviations ','interpreter','latex','Fontsize',15)
max_e3= max(max(e_list(:,3,:)))+0.01;
ylim([2,5.5 ])
xlim([dates(1),dates(end)]);
rshade(dates);
set(gcf,'paperpositionmode','manual','paperunits','inches');
dim = [10,8];
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure5_right');

