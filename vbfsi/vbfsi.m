function[Xiest,A,Sigma_A, B, Sigma_B_diag, J ,Sigma_J,wb,r]=vbfsi(Yf,Pf,A,Sigma_A,Bif,Sigma0,mu0,J,Sigma_J,rho)

%TUNED AND SCALED to 0.25
A_p=A;




Sigma_A_p=Sigma_A;
DIMRED_THR =0.5e5;
a_gamma0   = 1e-6;
b_gamma0   = 1e-6;
awb0= 1e-6; bwb0= 1e-6;
thr = 1e-4;
[m, n] = size(Yf);
obs=find(Pf==1);
p=length(obs)/(m*n);


[~,r]=size(Bif);
[~,S,~]=svd(A,'econ');

a_gamma=2*a_gamma0+m;

Xiest=A*Bif';

B=[mu0';Bif];
Y2sum = sum(Yf(:).^2);
scale2 = Y2sum / (m*n);
%scale = sqrt(scale2);
beta = 1./scale2;
%beta=1;

a_beta=p*m*n;
gammas = 1./diag(S.^2);
gammas = gammas(1:r);    %LOWER RAEE

sub=1;

wb=ones(1,r);

old_X=Xiest;

for it = 1:50
    
    % H=Y-E;
    H=Yf;
    
    W = diag(gammas);
    
    %% A step
    
    if sub==1
        
        JJ=J'*J+r*Sigma_J;
        psi_diag0=inv(Sigma0)+JJ;
        psi_diag=zeros(r,r,n);
        v1=zeros(n*r,1);
        %t=zeros(r,n);
        for j=1:n
            observed = find(Pf(:,j));
            Aj = A(observed,:);
            psi_diag(:,:,j)=eye(r)+JJ+beta*(Aj'*Aj+sum(Sigma_A(:,:,observed),3));
            v1((j-1)*r+1:j*r)=Aj'*H(observed,j);
            % t(:,j)=Aj'*H(observed,j);
            
        end
        % v1=reshape(t,n*r,1);
        psi_super_diag=-J';
        v0=(Sigma0)\mu0;
        
        % v1 =beta* reshape(A'*H,n*r,1);
        v=[v0; beta*v1];
        Sigma_B_diagT=zeros(r,r,n+1);
        Sigma_B_offdiagT=zeros(r,r,n);
        
        
        Sigma_B_diagT(:,:,1)=(psi_diag0)\eye(r);
        muBT=zeros((n+1)*r,1);
        muBT(1:r)=Sigma_B_diagT(:,:,1)*v(1:r);
        for i=1:n
            Sigma_B_offdiagT(:,:,i)=Sigma_B_diagT(:,:,i)*psi_super_diag;
            
            Sigma_B_diagT(:,:,i+1)=(psi_diag(:,:,i)-Sigma_B_offdiagT(:,:,i)'*psi_super_diag)\eye(r);
            
            muBT(i*r+1:(i+1)*r)=Sigma_B_diagT(:,:,i+1)*(v(i*r+1:(i+1)*r)-Sigma_B_offdiagT(:,:,i)'*muBT((i-1)*r+1:i*r));
            
        end;
        muB=zeros((n+1)*r,1);
        Sigma_B_diag(:,:,n+1)=Sigma_B_diagT(:,:,n+1);
        muB(n*r+1:(n+1)*r)=muBT(n*r+1:(n+1)*r);
        for i=n:-1:1
            Sigma_B_offdiag(:,:,i)=-Sigma_B_offdiagT(:,:,i)*Sigma_B_diag(:,:,i+1);
            Sigma_B_diag(:,:,i)=Sigma_B_diagT(:,:,i)-Sigma_B_offdiagT(:,:,i)*Sigma_B_offdiag(:,:,i)';
            muB((i-1)*r+1:i*r)=muBT((i-1)*r+1:i*r)-Sigma_B_offdiagT(:,:,i)*muB((i)*r+1:(i+1)*r);
            
        end;
        B=reshape(muB,[r n+1]); B=B';
        % end
        
        sub=2;
        Xiest = A*B(2:n+1,:)';
        
    else
        
        for i=1:m
            observed = find(Pf(i,:))+1;
            Bif = B(observed,:);
         
          % Sigma_A(:,:,i) = (beta*(Bif'*Bif) + beta*sum(Sigma_B_diag(:,:,observed),3) + W )\eye(r);
           % A(i,:) = (beta*(H(i,observed-1))*Bif+ A_p(i,:)*pinv(sum(Sigma_A(:,:,i),3)))*Sigma_A(:,:,i);
            
            
             Sigma_A(:,:,i) = (beta*(Bif'*Bif) + beta*sum(Sigma_B_diag(:,:,observed),3) + W+rho*pinv(sum(Sigma_A_p(:,:,i),3)+W.*1e-3 ) )\eye(r);
            A(i,:) = (beta*(H(i,observed-1))*Bif+ rho*A_p(i,:)*pinv(sum(Sigma_A_p(:,:,i),3)+W.*1e-3))*Sigma_A(:,:,i);
            
          % xyz=beta*(Bif'*Bif) + beta*sum(Sigma_A(:,:,i),3) + W;
           % Sigma_A(:,:,i) = (xyz)\eye(r);
            %A(i,:) = (beta*(H(i,observed-1))*Bif )*Sigma_A(:,:,i);
            %A(i,:) = (beta*(H(i,observed-1))*Bif + A(i,:)*pinv(Sigma_A(:,:,i)) )*Sigma_A(:,:,i);
            
        end
        Xiest = A*B(2:n+1,:)';
        sub=1;
        
    end
    
     
    b_beta=zeros(2,m);
    for l = 1: m
        observed = find(Pf(l, :))+1;
        b_beta(1,l) = sum(sum(B(observed,:).*(B(observed,:)*Sigma_A(:, :, l)))) ...
            + sum(sum(A(l, :).*(A(l, :) *sum(Sigma_B_diag(:, :, observed), 3)))) ...
            + sum(sum( Sigma_A(:, :, l).*sum(Sigma_B_diag(:, :, observed), 3)));
        
    end
    
     b_beta2=sum(b_beta(:))+sum(sum( abs(Yf - Pf.*(Xiest)).^2 ) );
    beta = (a_beta)/(b_beta2);

    if (sub==1  ) %CHECK-----------
        b_gamma = diag(A'*A)+ diag(sum(Sigma_A,3))+diag(B'*B) + diag(sum(Sigma_B_diag,3))+ b_gamma0;
        gammas = (a_gamma+n)./(b_gamma );
        MAX_GAMMA = min(gammas) * DIMRED_THR;
        
        if sum(find(gammas > MAX_GAMMA))
            indices = find(gammas <= MAX_GAMMA);
            
            A = A(:,indices);
            B = B(:,indices);
            A_p=A_p(:,indices);
            gammas = gammas(indices);
            
            Sigma_A = Sigma_A(indices,indices,:);
             Sigma_A_p = Sigma_A_p(indices,indices,:);
            %Sigma_B = Sigma_B(indices,indices);
            Sigma_B_diag=Sigma_B_diag(indices,indices,:);
            Sigma_B_offdiag=Sigma_B_offdiag(indices,indices,:);
            Sigma0=Sigma0(indices,indices);
            mu0=mu0(indices,:);
            wb=wb(indices);
            J=J(indices,indices);
            Sigma_J=Sigma_J(indices,indices);
            [m, r] = size(A);
        end
        
    end
    if sub==2
        Sigma_J=(diag(wb)+ B(1:n,:)'*B(1:n,:)+ sum(Sigma_B_diag(:,:,1:n),3))\eye(r);
        J=(Sigma_J*(B(1:n,:)'*B(2:n+1,:)+sum(Sigma_B_offdiag,3)))';
        wb=(2*awb0+ r)./(2*bwb0 + sum(J.^2)'+ r*diag(Sigma_J));
    end
    
    Xconv = sqrt(sum(sum(abs(old_X-Xiest).^2))/sum(sum(abs(old_X).^2)));
    fprintf('it %d: Xconv = %g, beta = %g, r = %d\n',it, Xconv,beta, r);
    % Check for convergence
    if it> 30 && Xconv < thr
        break;
    end
    old_X=Xiest;
    
end

end