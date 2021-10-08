clear all
clc
rng(1);
warning off;
%%
% dataset can be set to 1, 2,3,4
para.dataset =1;
% outlier percentage 0.05 and 0.1
    para.out_per=0.1;
    % error_v_1 is the error for vbfsi  and error_r_1 is the error rvbfsi
  [error_v_1,error_r_1] = Run_outlier_imputation(para);
  
