function [error_vbfsi,error_rvbfsi] = Run_outlier_imputation(para)
if isfield(para, 'dataset');      dataset = para.dataset;       else dataset = 1;   end
if isfield(para, 'out_per');      out_per = para.out_per;       else dataset = 0.1;   end




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
    rho=2;
     [m,n,d]=size(data);
    r=min(m,n);
    
end


%
fprintf("dataset is %d",dataset);
samp=[0.1,0.25,0.5,0.75];
[ss,ss2]=size(samp');
    error_vbfsi=0;
%%
for i=1:ss
    p=samp(i);
fprintf("sampling is %d",samp(i));

 fprintf("run RVBSF");
[mre_err,rmse_err]=rvbsf_run(data,p,start_day,end_day,rank,r,outlier_per);
[m1,m2]=size(mre_err);
     error_rvbsf(1:m2,i)=mre_err;
    error_rvbsf(1:m2,i+ss)=rmse_err;
    
 fprintf("run VBFSI");
    [mre_err,rmse_err]=vbfsi_run_outlier(data,p,start_day,end_day,rank,r,rho,out_per);
 [m1,m2]=size(mre_err);
    error_vbfsi(1:m2,i)=mre_err;
    error_vbfsi(1:m2,i+ss)=rmse_err;
    
fprintf("run RVBFSI");
     [mre_err,rmse_err]=rvbfsi_run(data,p,start_day,end_day,rank,r,rho,out_per);
 [m1,m2]=size(mre_err);
    error_rvbfsi(1:m2,i)=mre_err;
    error_rvbfsi(1:m2,i+ss)=rmse_err;
    


 
end


