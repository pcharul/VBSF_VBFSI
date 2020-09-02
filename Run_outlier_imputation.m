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
     load('/home/iiitd/AAAI_paper/datasets/air_data/air_quality_data.mat');
     start_day=30;
    end_day=60;
    rho=2;
     [m,n,d]=size(data);
    r=min(m,n);
    
end

error_trlf=0;
error_bcpf=0;
error_vmc=0;
error_vbsf=0;
error_1=0;
%
%
fprintf("dataset is %d",dataset);
samp=[0.05,0.1,0.15,0.25,0.5,0.75];
[ss,ss2]=size(samp');

%%
for i=1:ss
    p=samp(i);
fprintf("sampling is %d",samp(i));
 fprintf("run VBFSI");
    [mre_err,rmse_err]=vbfsi_run_outlier(data,p,start_day,end_day,rank,r,rho,out_per);
 [m1,m2]=size(mre_err);
    error_vbfsi(1:m2,i)=mre_err;
    error_vbfsi(1:m2,i+6)=rmse_err;
    
     [mre_err,rmse_err]=rvbfsi_run(data,p,start_day,end_day,rank,r,rho,out_per);
 [m1,m2]=size(mre_err);
    error_rvbfsi(1:m2,i)=mre_err;
    error_rvbfsi(1:m2,i+6)=rmse_err;


%  fprintf("run trlrf");
%    [mae_err,mape_err,rmse_err]=Run_TRLRF(data,p,start_day,end_day,rank,r);
% 
%     error_trlf(1:m2,i)=mae_err;
%     error_trlf(1:m2,i+6)=rmse_err;
% 
%     
%     %%
%      fprintf("run bcpf");
%     [mae_err,mape_err,rmse_err]=Run_BCPF(data,p,start_day,end_day,rank,r);  
% 
%     error_bcpf(1:m2,i)=mae_err;
%     error_bcpf(1:m2,i+6)=rmse_err;
% 
%     %%
%      p=samp(i);
% 
%  fprintf("run vmc");
%    [mae_err,mape_err,rmse_err]=Run_VMC(data,p,start_day,end_day,rank,r);
% 
%      error_vmc(1:m2,i)=mae_err;
%     error_vmc(1:m2,i+6)=rmse_err;
% 
%     %%
%      
%     %%
%      fprintf("run vbsf");
%    [mae_err,mape_err,rmse_err]=vbsf_run(data,p,start_day,end_day,rank,r);  
% 
%       error_vbsf(1:m2,i)=mae_err;
%     error_vbsf(1:m2,i+6)=rmse_err;
% 
%     
% 
%  fprintf("run trmf");
%    [mae_err,mape_err,rmse_err]=Run_TRMF(data,p,start_day,end_day,rank,r);
% 
%      error_trmf(1:m2,i)=mae_err;
% 
%     error_trmf(1:m2,i+6)=rmse_err;
    
    
end

%error=[
% csvwrite(strcat("result1/TRLRF1_error_data_",num2str(dataset),".csv"),error_trlf);
% csvwrite(strcat("result1/BCPF1_error_data_",num2str(dataset),".csv"),error_bcpf);
% 
% csvwrite(strcat("result1/VBSF1_error_data_",num2str(dataset),".csv"),error_vbsf);
% csvwrite(strcat("result1/VMC1_error_data_",num2str(dataset),".csv"),error_vmc);
% csvwrite(strcat("result1/TRMF1_error_data_",num2str(dataset),".csv"),error_1);
