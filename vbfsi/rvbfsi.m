function[Xiest,U,Sigma_U, V, Sigma_V_diag, F ,Sigma_F,wb,E]=rvbfsi(Yf,Pf,U,Sigma_U,Vif,Sigma0,mu0,F,Sigma_F,rho)
U_p=U;
Sigma_U_p=Sigma_U;


DIMRED_THR = 1e6;
thr = 1e-4;
[m, n] = size(Yf);
obs=find(Pf==1);
p=length(obs)/(m*n);


[~,r]=size(Vif);
[~,S,~]=svd(U,'econ');



Xiest=U*Vif';


V=[mu0';Vif];
Y2sum = sum(Yf(:).^2);
scale2 = Y2sum / (m*n);
scale = sqrt(scale2);
beta = 1./scale2;

a_beta=p*m*n;
gammas = 1./diag(S.^2);
gammas = gammas(1:r);   
Sigma_E = scale*ones(r,r);
alpha=ones(m,n)*scale;
E = randn(m, n) * sqrt(scale);
sub=1;
wb=ones(1,r);

old_X=Xiest;

for it = 1:50
   H=Yf-Pf.*E;
    W = diag(gammas);    
    if sub==1
        % Update V , F and respective hyperparameters 
        FF=F'*F+r*Sigma_F;
        psi_diag0=inv(Sigma0)+FF;
        psi_diag=zeros(r,r,n);
        v1=zeros(n*r,1);
        for j=1:n
            observed = find(Pf(:,j));
            Uj = U(observed,:);
            psi_diag(:,:,j)=eye(r)+FF+beta*(Uj'*Uj+sum(Sigma_U(:,:,observed),3));
            v1((j-1)*r+1:j*r)=Uj'*H(observed,j);
            
        end
        psi_super_diag=-F';
        v0=(Sigma0)\mu0;
        v=[v0; beta*v1];
        Sigma_V_diagT=zeros(r,r,n+1);
        Sigma_V_offdiagT=zeros(r,r,n);        
        Sigma_V_diagT(:,:,1)=(psi_diag0)\eye(r);
        muVT=zeros((n+1)*r,1);
        muVT(1:r)=Sigma_V_diagT(:,:,1)*v(1:r);
        for i=1:n
            Sigma_V_offdiagT(:,:,i)=Sigma_V_diagT(:,:,i)*psi_super_diag;            
            Sigma_V_diagT(:,:,i+1)=(psi_diag(:,:,i)-Sigma_V_offdiagT(:,:,i)'*psi_super_diag)\eye(r);            
            muVT(i*r+1:(i+1)*r)=Sigma_V_diagT(:,:,i+1)*(v(i*r+1:(i+1)*r)-Sigma_V_offdiagT(:,:,i)'*muVT((i-1)*r+1:i*r));            
        end
        muV=zeros((n+1)*r,1);
        Sigma_V_diag(:,:,n+1)=Sigma_V_diagT(:,:,n+1);
        muV(n*r+1:(n+1)*r)=muVT(n*r+1:(n+1)*r);
        for i=n:-1:1
            Sigma_V_offdiag(:,:,i)=-Sigma_V_offdiagT(:,:,i)*Sigma_V_diag(:,:,i+1);
            Sigma_V_diag(:,:,i)=Sigma_V_diagT(:,:,i)-Sigma_V_offdiagT(:,:,i)*Sigma_V_offdiag(:,:,i)';
            muV((i-1)*r+1:i*r)=muVT((i-1)*r+1:i*r)-Sigma_V_offdiagT(:,:,i)*muV((i)*r+1:(i+1)*r);
            
        end
        V=reshape(muV,[r n+1]); V=V';        
        sub=2;        
        Sigma_F=(diag(wb)+ V(1:n,:)'*V(1:n,:)+ sum(Sigma_V_diag(:,:,1:n),3))\eye(r);
        F=(Sigma_F*(V(1:n,:)'*V(2:n+1,:)+sum(Sigma_V_offdiag,3)))';
        wb=( r)./( sum(F.^2)'+ r*diag(Sigma_F));
        Xiest = U*V(2:n+1,:)';
        
    elseif sub==2
           % update U and gamma
            for i=1:m
                observed = find(Pf(i,:))+1;
                Vif = V(observed,:);
%                 psigma=pinv(sum(Sigma_U_p(:,:,i),3)+W );
%                 Sigma_U(:,:,i) = (beta*(Vif'*Vif) + beta*sum(Sigma_V_diag(:,:,observed),3) +W + rho*psigma )\eye(r);
%                 U(i,:) = (beta*(H(i,observed-1))*Vif+ rho*U_p(i,:)*pinv(sum(Sigma_U_p(:,:,i),3)))*Sigma_U(:,:,i);
            Sigma_U(:,:,i) = (beta*(Vif'*Vif) + beta*sum(Sigma_V_diag(:,:,observed),3) + W+rho*pinv(sum(Sigma_U_p(:,:,i),3)+1e-4) )\eye(r);
            U(i,:) = (beta*(H(i,observed-1))*Vif+ rho*U_p(i,:)*pinv(sum(Sigma_U_p(:,:,i),3)+1e-4))*Sigma_U(:,:,i);
            

            end
            Xiest = U*V(2:n+1,:)';
       
        sub=3;
        % update gamma and rank determination
        
        b_gamma =diag(V'*V) + diag(sum(Sigma_V_diag,3)) + diag(U'*U)+ diag(sum(Sigma_U,3));
        gammas = (m + n)./(b_gamma );
        MAX_GAMMA = min(gammas) * DIMRED_THR/2;
        
        if sum(find(gammas > MAX_GAMMA))
            indices = find(gammas <= MAX_GAMMA);            
            U = U(:,indices);
            V = V(:,indices);
            U_p=U_p(:,indices);
            gammas = gammas(indices);            
            Sigma_U = Sigma_U(indices,indices,:);
            Sigma_U_p = Sigma_U_p(indices,indices,:);   
            Sigma_V_diag=Sigma_V_diag(indices,indices,:);
            Sigma_V_offdiag=Sigma_V_offdiag(indices,indices,:);
            Sigma0=Sigma0(indices,indices);
            mu0=mu0(indices,:);
            wb=wb(indices);
            F=F(indices,indices);
            Sigma_F=Sigma_F(indices,indices);
            [m, r] = size(U);
            
        
        end
        
    else
        % update E and alpha
          Sigma_E=1./(alpha +beta);
        E=beta*(Pf.*Yf -Pf.*(Xiest)).*Sigma_E;
        alpha= (1-alpha.*Sigma_E )./ (E.^2 + eps );
        sub=1;
    end
    b_beta=zeros(2,m);
    for l = 1: m
        observed = find(Pf(l, :))+1;
        b_beta(1,l) = sum(sum(V(observed,:).*(V(observed,:)*Sigma_U(:, :, l)))) ...
            + sum(sum(U(l, :).*(U(l, :) *sum(Sigma_V_diag(:, :, observed), 3)))) ...
            + sum(sum( Sigma_U(:, :, l).*sum(Sigma_V_diag(:, :, observed), 3)));
        
    end
    
     b_beta2=sum(b_beta(:))+sum(sum( abs(Yf - Pf.*(Xiest)-Pf.*E).^2 ) )+sum(sum((Sigma_E)));
    beta = (a_beta)/(b_beta2);
    Xconv = sqrt(sum(sum(abs(old_X-Xiest).^2))/sum(sum(abs(old_X).^2)));
    if it> 50 && Xconv < thr
        fprintf('break')
        break;
    end
    old_X=Xiest;
    
end

end
