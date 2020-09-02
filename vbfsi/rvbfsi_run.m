function [mre_err,rmse_err] = rvbfsi_run(data,p,start_day,end_day,rank,r,rho,s)
A_1=0;
P_1=0;

if rho==1

rho=1.146*exp(-4.1*p)+0.0438*exp(0.8*p);
elseif rho ==2
    
rho=1.282*exp(-11.18*p)+0.0289*exp(1.74*p);
else
    rho=0.2;
end
for iij=start_day-8:start_day-1
    Ycheck=data(:,:,iij);
    [m, n]=size(Ycheck);
   
    Obs1 = randperm(m*n); Obs = Obs1(1:round(p*m*n));
    P = zeros(m,n);  P(Obs) = 1;
   slot=n;
    Y= P.*Ycheck;
    Y=add_outlier_mtx(Y,s);
if rank>100
   r=fix(min(m,n)/2);
end
    Yi=Y;
      
        Pi=P;
        
        [U, S, V] = svd(Yi, 'econ');
        A = U(:,1:r)*(S(1:r,1:r).^(0.5));
        Y2sum = sum(Yi(:).^2);
        scale2 = Y2sum / (m*slot);
        scale = sqrt(scale2);
        [m, r]=size(A);
        sig_A=scale *repmat(eye(r),[1 1 m]);
        Sigma_A=sig_A;
        B = (S(1:r,1:r)).^(0.5)*V(:,1:r)';
        B = B';
        sigma0=1e3*eye(r);
        Sigma0=sigma0;
        mu0=zeros(r,1);


        J=zeros(r,r); Sigma_J=eye(r,r);
        [Xiest,A,Sigma_A, B, Sigma_B_diag, J ,Sigma_J,wb,r]=rvbfsi(Yi,Pi,A,Sigma_A,B,Sigma0,mu0,J,Sigma_J,0);
     
      A_1=A_1+A;
    
end
A=A_1./8;

%%
for iij=start_day:end_day
  
     Yi_true=data(:,:,iij);
    [m, n]=size(Yi_true);
  
   
    Obs1 = randperm(m*n); Obs = Obs1(1:round(p*m*n));
 
    P = zeros(m,n);  P(Obs) = 1;
   slot=n;
 
  

    Yi=P.*Yi_true;
       Yi=add_outlier_mtx(Yi,s);
       
        Pi=P;
        
        [U, S, V] = svd(Yi, 'econ');
      
        Y2sum = sum(Yi(:).^2);
        scale2 = Y2sum / (m*slot);
        scale = sqrt(scale2);
        [m, r]=size(A);

        B = (S(1:r,1:r)).^(0.5)*V(:,1:r)';
        B = B';
        sigma0=1e3*eye(r);
        Sigma0=sigma0;
        mu0=zeros(r,1);    
         [Xiest,A,Sigma_A, B_1, Sigma_B_diag_1, J ,Sigma_J,wb_1,r]=rvbfsi(Yi,Pi,A,Sigma_A,B,Sigma0,mu0,J,Sigma_J,rho);
        
    mre_err(iij-start_day+1)=mre_error(Yi_true,Xiest,P);  
   
    rmse_err(iij-start_day+1)=rmse_error(Yi_true,Xiest,P);  
        fprintf("day is %d rmse is is %d\n",iij,rmse_error(Yi_true,Xiest,P));        
end
end

