clear all
clc
rng(1);
para.dataset =1;
    para.out_per=0.05;
  [error_v_1,error_r_1] = Run_outlier_imputation(para);
  para.dataset =1;

    para.out_per=0.1; 
  [error_v_2,error_r_2] = Run_outlier_imputation(para);
  para.dataset =1;
        para.out_per=0.02; 
   [error_v_3,error_r_3]  = Run_outlier_imputation(para);
  