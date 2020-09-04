clear all
clc
rng(1);
para.dataset =1;
    
  error_1 = Run_imputation(para);
  para.dataset =2;
    
  error_2 = Run_imputation(para);
  para.dataset =3;
    
  error_3 = Run_imputation(para);
  csvwrite("error_1.csv",error_1);
   csvwrite("error_2.csv",error_2);
    csvwrite("error_3.csv",error_3);