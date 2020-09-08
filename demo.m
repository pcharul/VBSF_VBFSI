clear all
 clc
 warning off;
rng(1);
% para.dataset =1;
%     
%   error_1 = Run_imputation(para);
   para.dataset =1;
%     
  error_vbfsi = Run_imputation(para);
%   para.dataset =3;
%     
%   error_3 = Run_imputation(para);
   csvwrite("error_act_vbfsi.csv",error_vbfsi);
%    csvwrite("error_2.csv",error_2);
%     csvwrite("error_3.csv",error_3);





    para.out_per=0.1;
  [error_v_1,error_r_1] = Run_outlier_imputation(para);
%   para.dataset =1;
% 
    para.out_per=0.05; 
  [error_v_2,error_r_2] = Run_outlier_imputation(para);
%   para.dataset =1;
%         para.out_per=0.02; 
%    [error_v_3,error_r_3]  = Run_outlier_imputation(para);
%   
%    
   %csvwrite(strcat("error_v_1.csv"),error_v_1,1);
   csvwrite(strcat("error_r_5_1.csv"),error_r_2);
   csvwrite(strcat("error_v_5_1.csv"),error_v_2);
   
   csvwrite(strcat("error_r_1_1.csv"),error_r_1);
   csvwrite(strcat("error_v_1_1.csv"),error_v_1);