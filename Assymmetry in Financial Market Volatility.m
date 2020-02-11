clear all
close all
clc

load('crsp_vwret.mat');

[yr,mo,day,dates] = decode_date(crspdates);
Theta0_gjr = [0 0.003 0.3 0.6 0.1]; 

A_gjr = [0 0 1 1 0.5];
b_gjr = 1;
LB_gjr = [-1 0 0 0 0];
UB_gjr = [];

Theta_gjr = fmincon(@(Theta) gjrgarchlgl(Theta,ret),Theta0_gjr,A_gjr,b_gjr,[],[],LB_gjr,UB_gjr,[]); % gjrgarch1,1

gjr_lgl = gjrgarchlgl(Theta_gjr,ret);

S = score_matrix(@gjrgarchlgl,Theta_gjr,ret); 
H = hessian_2sided(@gjrgarchlgl,Theta_gjr,ret);

vcov = H\S/H;
seQMLE = diag(sqrt(vcov))';

tstat = Theta_gjr./seQMLE;

disp('Estimation results')
Theta_gjr
seQMLE
tstat

% Constraints
% 1. Upper and lower bound of parameters
LB = [-1 0 0 0];
UB = [];
% 2. Sum constraints (alpha + beta) < 1
A = [0 0 1 1];
b = 1;

options = optimset('fmincon');
options.Display = 'none';
warning off

% Use estimated parameters as starting values
Theta0 = [0 0.003 0.3 0.6]; 
 
Theta = fmincon(@(Theta) restrictedGARCHlgl(Theta,ret),Theta0, A,b,[],[],LB,UB,[]);

garch_lgl = restrictedGARCHlgl(Theta, ret);

LRT = 2*((-gjr_lgl) - (-garch_lgl))

chi2_crit = chi2inv(1-0.05,1)

LRT > chi2_crit; 
[lgl_gjr, lglt_gjr, z_gjr] = gjrgarchlgl(Theta_gjr,ret);
sig_gjr = sqrt(z_gjr);

[lgl_garch, lglt_garch, z_garch] =  restrictedGARCHlgl(Theta, ret);
sig_garch = sqrt(z_garch);

figure;
hold on
plot(dates,abs(ret),'y--','LineWidth',2), datetick
plot(dates,sig_gjr,'r',dates,sig_garch,'b','LineWidth',2), datetick
legend('Returns','GJR','GARCH')
hold off
figure;
plot(dates,sig_gjr,'r',dates,sig_garch,'b'), datetick
legend('GJR','GARCH')
