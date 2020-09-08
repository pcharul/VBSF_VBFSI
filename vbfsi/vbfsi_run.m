function [mre_err,rmse_err] = vbfsi_run(data,p,start_day,end_day,rank,r,rho)
Y_1=0;
P_1=0;



%rho=1.146*exp(-4.1*p)+0.0438*exp(0.8*p);
%rho=1.167*exp(-4.334*p);
%rho=1.29 *exp(-4.806*p) + 0.0002809*exp(7.289*p);
% elseif rho ==2
%     
 %rho=1.282*exp(-11.18*p)+0.0289*exp(1.74*p);
% else
%     rho=0.2;
% end

%rho=-0.35714*p + 0.26786;


%rho=1.094*exp(-3.871*p) + 0.008622*exp(3.764*p)-0.1;


%rho=1.146*exp(-4.1*p)+0.0438*exp(0.8*p);
%rho=1.167*exp(-4.334*p);
%rho=1.29 *exp(-4.806*p) + 0.0002809*exp(7.289*p);
if rho==1
    rho=1.094*exp(-3.871*p) + 0.008622*exp(3.764*p);
 elseif rho ==2
%     
 rho=1.282*exp(-11.18*p)+0.0289*exp(1.74*p);
 else
     rho=0.2;
 end

%rho=-0.35714*p + 0.26786;


%rho=1.094*exp(-3.871*p) + 0.008622*exp(3.764*p)-0.1;


for iij=start_day-8:start_day-1
    Ycheck=data(:,:,iij);
    [m, n]=size(Ycheck);
   
   

  
    Obs1 = randperm(m*n); Obs = Obs1(1:round(p*m*n));
    P = zeros(m,n);  P(Obs) = 1;
   slot=n;
 
    Y= P.*Ycheck;
    Y_1=Y_1+Y;
    P_1=P_1+P;
end
Y=Y_1./P_1;
Y(isnan(Y))=0;
%%
  P=P_1./P_1;
  
P(isnan(P))=0;
%% 
  
  
if r>100
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
        [Xiest,A,Sigma_A, B, Sigma_B_diag, J ,Sigma_J,wb,r]=vbfsi(Yi,Pi,A,Sigma_A,B,Sigma0,mu0,J,Sigma_J,0);
     
      
 

for iij=start_day:end_day
  
     Yi_true=data(:,:,iij);
    [m, n]=size(Yi_true);
  
   %   Obs=RT(:,:,iij);
  %  P = round(Obs+0.5-p);
  
    Obs1 = randperm(m*n); Obs = Obs1(1:round(p*m*n));
    P = zeros(m,n);  P(Obs) = 1;
   slot=n;
 


    Yi=P.*Yi_true;
        
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
         [Xiest,A,Sigma_A, B_1, Sigma_B_diag_1, J ,Sigma_J,wb_1,r]=vbfsi(Yi,Pi,A,Sigma_A,B,Sigma0,mu0,J,Sigma_J,rho);
        
    mre_err(iij-start_day+1)=mre_error(Yi_true,Xiest,P);  
 
    rmse_err(iij-start_day+1)=rmse_error(Yi_true,Xiest,P);
end



