clear all
close all 
clc
%%
[Data, text, alldata] = xlsread('DataH2MSF6321.xlsx');

Ret=Data(:,2:7); % These are the returns of the 6 assets
N=size(Ret,2);
E0=mean(Ret,1)'; % Here E0 is just the in sample mean of the EXCESS returns
SIGMA=cov(Ret);

RF=0.4;%

E0=E0+RF; % So that we move from excess returns of returns

% Get the predictors

DP=Data(:,8);
TERM=Data(:,9);
DEF=Data(:,10);
INF=Data(:,11);
UNE=Data(:,12);


  

%%
% Just for fun, here are some basic characteristics of the assets
SDs=sqrt(diag(SIGMA))';
SRs=(E0'-RF)./SDs;

Table_Characteristics=[E0';SDs;SRs];
disp(' '); 
disp(' Mean, Sdev and Sharpe ratio of the 6 assets')
row= {'Mean','Sdev','Sharpe ratio'};
head= {'Asset 1','Asset 2','Asset 3','Asset 4','Asset 5','Asset 6'};
fmt = ' %6.2f ';
texprint(Table_Characteristics, fmt, head, row);

% ***********************************************************************
% Part a) Run a time series regression of returns on predictors. Important:
% the predictor needs to be lagged one period
% ***********************************************************************

%% Question a - OLS Estimates 

T= 6;
summarytable=[];
rsquares=[];
for i = 1:T
    table= fitlm(DP(1:end-1),Ret(2:end,i));
    temp=table2array(table.Coefficients);
    summarytable=[summarytable;temp];
    rsquares=[rsquares; table.Rsquared.Ordinary];
end
alphas = summarytable([1 3 5 7 9 11],1)
betas = summarytable([2 4 6 8 10 12],1)
tstats = summarytable(:,3)

% compute expected return of 6 assets for 2006-07
E=[];
for i = 1:T
    temp = RF + (alphas(i)+ betas(i)*DP(end));
    E = [E; temp];
end

%% Question B - MVEP without Short Sell

% Using optimizers (very useful when you want to calculate optimal portfolios with constraints)

fun=@(w)-(w'*E-RF)/sqrt(w'*SIGMA*w); % note I introduce a minus here because we want to maximize, not minimize (which is what the fmincon function does)

%options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options = optimoptions('fmincon','OptimalityTolerance',1e-6,'StepTolerance',1e-6,'FunctionTolerance', 1e-6,'ConstraintTolerance',1e-6,'MaxFunctionEvaluations', 2000,'HonorBounds',true);
% the Mean Variance Efficient Portfolio (that is, the portfolio with Maximum Sharpe ratio)
    w0=1/6*ones(6,1);
%  No short sale constraints
    [w fvalSR0] = fmincon(fun,w0,[],[],[ones(1,N)],[1])
    
    ERp = w'*E;
    SDp = sqrt(w'*SIGMA*w);
 
    w
    % W[MVEP]: optimal weights for MVEP 
    % ER[MVEP]: Wi * Mean(Ret(i))     
    % SD[MVEP]: Standard Deviation for MVEP

    % Observation: 
    % Weights are not reasonable as you have negative weights, suggesting
    % that investors are short-selling their investments 
%% Question C - MVEP With short sale constraints
    [w_shortsell fvalSR1] = fmincon(fun,w0,[],[],[ones(1,N)],[1],zeros(N,1), ones(N,1),[],options)

    SRp = [fvalSR1./fvalSR0 fvalSR1 fvalSR0]'
    % Observation:
    % see an increase in SR from imposing a short sale constraint on the
    % optimization model 

%% Question D - Redo A,B,C w/ more predictors 

summarytable2=[];
rsquares2=[];
X= [DP(1:end-1) TERM(1:end-1) DEF(1:end-1) INF(1:end-1) UNE(1:end-1) ];
for i = 1:T
    table= fitlm(X,Ret(2:end,i));
    temp=table2array(table.Coefficients);
    summarytable2=[summarytable2;temp];
    rsquares2=[rsquares2; table.Rsquared.Ordinary];
end
 coef = summarytable2(:,1);
 tstats2 = summarytable2(:,3);
 % R-square with more predictors seem to produce a better fit in
 % prediciting returns than just D/P ratio
% compute expected return of 6 assets for 2006-07
E2=[];

for i =1:6:length(coef)
    temp = RF + (coef(i)+ coef(i+1)*DP(end)+ coef(i+2)*TERM(end)+ coef(i+3)*DEF(end) +coef(i+4)*INF(end) +coef(i+5)*UNE(end));
    E2 = [E2; temp];
end

fun=@(w)-(w'*E2-RF)/sqrt(w'*SIGMA*w); % note I introduce a minus here because we want to maximize, not minimize (which is what the fmincon function does)

%options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options = optimoptions('fmincon','OptimalityTolerance',1e-6,'StepTolerance',1e-6,'FunctionTolerance', 1e-6,'ConstraintTolerance',1e-6,'MaxFunctionEvaluations', 2000,'HonorBounds',true);
% the Mean Variance Efficient Portfolio (that is, the portfolio with Maximum Sharpe ratio)
    w0=1/6*ones(6,1);
%  No short sale constraints                      %lower bond %upper bond
[w2 fval] = fmincon(fun, w0,[],[],[ones(1,N)],[1],zeros(6,1),ones(6,1),[],options);
    
    ERp2 = w2'*E2;
    SDp2 = sqrt(w2'*SIGMA*w2);

% The MVEP mean return decreases as well as the standard errors 


%% Question E - Create Minimum Variance Frontier 
% Set Up 
EM = mean(Ret)';
A=EM'*inv(SIGMA)*EM;
B=EM'*inv(SIGMA)*ones(N,1);

C=ones(1,N)*inv(SIGMA)*ones(N,1);
w_minvar=inv(SIGMA)*ones(N,1)/C;
mu_minvar=B/C;
SD_minvar=sqrt(1/C);

mu=[-5:0.05:5]; % that is, I will consider a target return from -10% to +10% 

MEAN_SD1=[];
for i=1:length(mu) % this loop computes the minimum variance (standard deviation) for each target return
    var_p=(C*mu(i)^2-2*B*mu(i)+A)/(A*C-B^2);
    tmp=[mu(i),sqrt(var_p)];
    MEAN_SD1=[MEAN_SD1;tmp];
end

% i) sample mean 
hold on
p1=plot(MEAN_SD1(:,2),MEAN_SD1(:,1),'color','red','LineWidth',1.5); % minimum-variance frontier
xlabel({'Standard Deviation (in %)'});
ylabel({'Expected Return (in %)'});
hold off

% Set Up 
A=E'*inv(SIGMA)*E;
B=E'*inv(SIGMA)*ones(N,1);

C=ones(1,N)*inv(SIGMA)*ones(N,1);
w_minvar=inv(SIGMA)*ones(N,1)/C;
mu_minvar=B/C;
SD_minvar=sqrt(1/C);

mu=[-5:0.05:5]; % that is, I will consider a target return from -10% to +10% 

MEAN_SD2=[];
for i=1:length(mu) % this loop computes the minimum variance (standard deviation) for each target return
    var_p=(C*mu(i)^2-2*B*mu(i)+A)/(A*C-B^2);
    tmp=[mu(i),sqrt(var_p)];
    MEAN_SD2=[MEAN_SD2;tmp];
end

% ii) given by the predicted values using equation (1)
hold on
p2=plot(MEAN_SD2(:,2),MEAN_SD2(:,1),'color','black','LineWidth',1.5); % minimum-variance frontier
xlabel({'Standard Deviation (in %)'});
ylabel({'Expected Return (in %)'});
hold off

% Set Up 
A=E2'*inv(SIGMA)*E2;
B=E2'*inv(SIGMA)*ones(N,1);

C=ones(1,N)*inv(SIGMA)*ones(N,1);
w_minvar=inv(SIGMA)*ones(N,1)/C;
mu_minvar=B/C;
SD_minvar=sqrt(1/C);

mu=[-5:0.05:5]; % that is, I will consider a target return from -10% to +10% 

MEAN_SD3=[];
for i=1:length(mu) % this loop computes the minimum variance (standard deviation) for each target return
    var_p=(C*mu(i)^2-2*B*mu(i)+A)/(A*C-B^2);
    tmp=[mu(i),sqrt(var_p)];
    MEAN_SD3=[MEAN_SD3;tmp];
end

% iii) given by the predicted values given by equation (2)
hold on
p3= plot(MEAN_SD3(:,2),MEAN_SD3(:,1),'color','blue','LineWidth',1.5); % minimum-variance frontier
xlabel({'Standard Deviation (in %)'});
ylabel({'Expected Return (in %)'});
hold off

legend([p1 p2 p3],'Red: Sample Mean', 'Black: Predicted Values using Equation (1)', 'Blue: Predicted Values using Equation (2)','Location', 'east')
title('Minimum Variance Frontier Plot ')
