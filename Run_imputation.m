
function [error_vbfsi] = Run_imputation(para)
if isfield(para, 'dataset');      dataset = para.dataset;       else dataset = 1;   end

if isfield(para, 'rho');      rho = para.rho;       else rho = 3;   end

rank=1;
r=0;
if dataset==1
    load('datasets/delhi.mat');
    start_day=30;
    end_day=60;
    rho=1;
    [m,n,d]=size(data);
    r=min(m,n);
elseif dataset==2
    load('datasets/pems_data.mat');
    start_day=16;
    end_day=44;
    rho=1;
    
     [m,n,d]=size(data);
    r=min(m,n);
elseif dataset==3
    load('datasets/tensor.mat');
start_day=30;
    end_day=60;
    rho=1;
     [m,n,d]=size(data);
    r=min(m,n);
elseif dataset==4
     load('datasets/air_quality_data.mat');
     start_day=30;
    end_day=60;
   
     [m,n,d]=size(data);
    r=min(m,n);
    
end

fprintf("dataset is %d",dataset);
samp=[0.05,0.15,0.25,0.5,0.75];

[ss,ss2]=size(samp');

%%
for i=1:5
    p=samp(i);

fprintf("sampling is %d",samp(i));
 fprintf("run VBFSI");
 
  fprintf("run vbsf");
 [mre_err,rmse_err]=vbsf_run(data,p,start_day,end_day,rank,r);  
   [m1,m2]=size(mre_err);
      error_vbsf(1:m2,i)=mre_err;
    error_vbsf(1:m2,i+5)=rmse_err;

    [mre_err,rmse_err]=vbfsi_run(data,p,start_day,end_day,rank,r,rho);
 [m1,m2]=size(mre_err);
    error_vbfsi(1:m2,i)=mre_err;
    error_vbfsi(1:m2,i+5)=rmse_err;

    
 

    
end

