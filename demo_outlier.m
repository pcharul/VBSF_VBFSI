clear all
clc
rng(1);
warning off;
para.dataset =1;
    para.out_per=0.1;
  [error_v_1,error_r_1] = Run_outlier_imputation(para);
%   para.dataset =1;
% 
%     para.out_per=0.1; 
%   [error_v_2,error_r_2] = Run_outlier_imputation(para);
%   para.dataset =1;
%         para.out_per=0.02; 
%    [error_v_3,error_r_3]  = Run_outlier_imputation(para);
%   
%    
   %csvwrite(strcat("error_v_1.csv"),error_v_1,1);
   csvwrite(strcat("error_r_1.csv"),error_r_1);
%    csvwrite(strcat("error_v_2.csv"),mean(error_v_2,1));
%    csvwrite(strcat("error_r_2.csv"),mean(error_r_2,1));
%    csvwrite(strcat("error_v_3.csv"),mean(error_v_3,1));
%    csvwrite(strcat("error_r_3.csv"),mean(error_r_3,1));