function[Xiest,A,Sigma_A, B, Sigma_B_diag, J ,Sigma_J,wb,E]=rvbfsi(Yf,Pf,A,Sigma_A,Bif,Sigma0,mu0,J,Sigma_J,rho)
A_p=A;
Sigma_A_p=Sigma_A;
%TUNED AND SCALED to 0.25

DIMRED_THR = 1e6;
a_gamma0   = 0;
b_gamma0   = 0;
awb0= 0; bwb0=0;
thr = 1e-4;
[m, n] = size(Yf);
obs=find(Pf==1);
p=length(obs)/(m*n);


[~,r]=size(Bif);
[~,S,~]=svd(A,'econ');



Xiest=A*Bif';


B=[mu0';Bif];
Y2sum = sum(Yf(:).^2);
scale2 = Y2sum / (m*n);
scale = sqrt(scale2);
beta = 1./scale2;
%beta=1;

a_beta=p*m*n;
gammas = 1./diag(S.^2);
gammas = gammas(1:r);    %LOWER RAEE
Sigma_E = scale*ones(r,r);
alpha=ones(m,n)*scale;
%E=(Yf - Pf.*(Xiest));
E = randn(m, n) * sqrt(scale);
sub=1;
wb=ones(1,r);

old_X=Xiest;

for it = 1:50
    
    % H=Y-E;
   H=Yf-Pf.*E;
    %fprintf('iteration is %d',it);
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
        
        Sigma_J=(diag(wb)+ B(1:n,:)'*B(1:n,:)+ sum(Sigma_B_diag(:,:,1:n),3))\eye(r);
        J=(Sigma_J*(B(1:n,:)'*B(2:n+1,:)+sum(Sigma_B_offdiag,3)))';
        wb=(2*awb0+ r)./(2*bwb0 + sum(J.^2)'+ r*diag(Sigma_J));
        Xiest = A*B(2:n+1,:)';
        
    elseif sub==2
        
        for i=1:m
            observed = find(Pf(i,:))+1;
            Bif = B(observed,:);
            psigma=pinv(sum(Sigma_A_p(:,:,i),3)+W );
            Sigma_A(:,:,i) = (beta*(Bif'*Bif) + beta*sum(Sigma_B_diag(:,:,observed),3) +W + rho*psigma )\eye(r);
             %Sigma_A(:,:,i) = (beta*(Bif'*Bif) + beta*sum(Sigma_B_diag(:,:,observed),3) + W)\eye(r);
            A(i,:) = (beta*(H(i,observed-1))*Bif+ rho*A_p(i,:)*pinv(sum(Sigma_A_p(:,:,i),3)))*Sigma_A(:,:,i);
            
            
        end
        Xiest = A*B(2:n+1,:)';
      
        
  
      
        
        sub=3;
      
        
       b_gamma =diag(B'*B) + diag(sum(Sigma_B_diag,3)) + diag(A'*A)+ diag(sum(Sigma_A,3))+ b_gamma0;
        gammas = (m + n+ a_gamma0)./(b_gamma );
        MAX_GAMMA = min(gammas) * DIMRED_THR/2;
        
        if sum(find(gammas > MAX_GAMMA))
            %fprintf('gamma')
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
        
    else
          Sigma_E=1./(alpha +beta);
        E=beta*(Pf.*Yf -Pf.*(Xiest)).*Sigma_E;

        %alpha = 1./(E.^2 + Sigma_E);
        alpha= (1-alpha.*Sigma_E )./ (E.^2 + eps );
        sub=1;
        
    %err = err + trace(B(observed,:)'*B(observed,:)*Sigma_A(:, :, l)) ...
                %+ trace(A(l, :)'*A(l, :) *sum(Sigma_B(:, :, observed), 3)) ...
                %+ trace( Sigma_A(:, :, l)*sum(Sigma_B(:, :, observed), 3));
                
                
     %err = sum(sum( abs(Y - X - E).^2 ) ) + n*trace(A'*A*Sigma_B) + m*trace(B'*B*Sigma_A) + m*n*trace(Sigma_A*Sigma_B) + sum(sum((Sigma_E)));
     
    b_beta=zeros(2,m);
    for l = 1: m
        observed = find(Pf(l, :))+1;
        b_beta(1,l) = sum(sum(B(observed,:).*(B(observed,:)*Sigma_A(:, :, l)))) ...
            + sum(sum(A(l, :).*(A(l, :) *sum(Sigma_B_diag(:, :, observed), 3)))) ...
            + sum(sum( Sigma_A(:, :, l).*sum(Sigma_B_diag(:, :, observed), 3)));
        
    end
    
     b_beta2=sum(b_beta(:))+sum(sum( abs(Yf - Pf.*(Xiest)-Pf.*E).^2 ) )+sum(sum((Sigma_E)));
    beta = (a_beta)/(b_beta2);
   

        
    end


    Xconv = sqrt(sum(sum(abs(old_X-Xiest).^2))/sum(sum(abs(old_X).^2)));
    if it> 50 && Xconv < thr
        fprintf('break')
        break;
    end
    old_X=Xiest;
    
end

end
