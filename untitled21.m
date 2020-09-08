clear all
clc
rng(1);

rank=1;
r=0;
 load('datasets/air_quality_data.mat');
     start_day=30;
    end_day=60;
    rho=2;
     [m,n,d]=size(data);
    r=min(m,n);
 samp=[0.05,0.1,0.15,0.25,0.5,0.75];

[ss,ss2]=size(samp');
rho=0;
%%
for i=1:7
    p=samp;
  
fprintf("sampling is %d",samp(i));
 fprintf("run VBFSI");
    [mre_err,rmse_err]=vbfsi_run(data,p,start_day,end_day,rank,r,rho);
 [m1,m2]=size(mre_err);
    error_vbfsi(1,i)=mean(mre_err);
    error_vbfsi(2,i)=mean(rmse_err);

end