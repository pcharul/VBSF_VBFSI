clear all
 clc
 addpath(genpath(pwd));
 warning off;
 rng('default')
rng(1);
% dataset can be set to 1, 2,3,4
para.dataset =1;
% % rho=1 for traffic data i.e dataset 1,2,3 and rho=2 for air quality data i.e. dataset=4.     
para.rho=2;
error_1 = Run_imputation(para);
error_1_m=mean(error_1,1);

