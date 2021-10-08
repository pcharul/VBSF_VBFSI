function [mre_err,rmse_err] = rvbsf_run(data,p,start_day,end_day,rank,r,s)

for iij=start_day: end_day
    fprintf('day no %d: \n',iij);
     Yi_true=data(:,:,iij);
    [m, n]=size(Yi_true);

  Yi_true(isnan(Yi_true))=0;
    Obs1 = randperm(m*n); Obs = Obs1(1:round(p*m*n));
    P = zeros(m,n);  P(Obs) = 1;
   slot=n;

    Y= P.*Yi_true;
    Y=add_outlier_mtx(Y,s);
if rank==1
    r1=r;
else
   r1=min(m,n);
end
r=r1;

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

       
     
        
         [Xiest,A,Sigma_A, B_1, Sigma_B_diag_1, J ,Sigma_J,wb_1,r]=rvbsf(Yi,Pi,A,Sigma_A,B,Sigma0,mu0,J,Sigma_J);
  mre_err(iij-start_day+1)=mre_error(Yi_true,Xiest,P);     
    rmse_err(iij-start_day+1)=rmse_error(Yi_true,Xiest,P);  
         
 
end
end

