function ETI = Cal_ETI(g_mu,g_sig,Gbest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the ETI(x|w) at multiple requests, e.g., n  
% g_mu  : n*M  predictive mean
% g_sig : n*M  square root of the predictive variance
% Gbest : n*1   
%------------------------------- Reference --------------------------------
% L. Zhao and Q. Zhang, Exact Formulas  for the Computation of Expected 
% Tchebycheff Improvement. Proceedings of the IEEE Congress on Evolutionary 
% Computation,2023.
%--------------------------------------------------------------------------
% This function is written by Liang Zhao
% https://github.com/mobo-d/MOEAD-EGO
    
    [n,M] = size(g_mu); ETI = zeros(n,1);
    g_sig(g_sig<0) = 0; g_sig2 = g_sig.^2; % n*M
    
    if M == 2
        tau = (Gbest-g_mu)./g_sig; 
        Cdf = normcdf(tau);  Pdf = normpdf(tau); 
        sig_m = sqrt(sum(g_sig2,2));
        diff_mu = (g_mu(:,1)-g_mu(:,2));
        tau12 = diff_mu./sig_m;
        norm_sig = g_sig./sig_m;

        ETI = g_sig(:,1).*Pdf(:,1).*Cdf(:,2)+g_sig(:,2).*Pdf(:,2).*Cdf(:,1)+(Gbest-g_mu(:,1)).* prod(Cdf,2)-...
            sig_m.*normpdf(tau12).*normcdf(tau(:,2).*norm_sig(:,1)+tau(:,1).*norm_sig(:,2));

        for i = 1:1:n
            BVN_term = bvnl(-tau12(i),tau(i,2),-norm_sig(i,2));
            ETI(i) = ETI(i) + diff_mu(i).*BVN_term; 
        end
    elseif M > 2
        for k=1:M
            L_k =  eye(M);    L_k(:,k)=-1;  L_k(k,k)= 1;
            mu_k = g_mu*L_k'; % The mean of Z^(k), n*m
            b_k = zeros(n,M); b_k(:,k) = Gbest;
            for j=1:1:n 
                  Sigma_k =L_k*diag(g_sig2(j,:))*L_k'; % The covariance matrix of Z^(k) for j-th query, m*m
                  % P(Z^k<=b^k)
                  p_k  =   qsimvnv(200*M,Sigma_k, -inf.*ones(M,1), (b_k(j,:)-mu_k(j,:))') ; 
                  ETI(j) =  ETI(j)+ p_k *(Gbest(j)-mu_k(j,k));

                  temp_ci = (b_k(j,:) -mu_k(j,:))'./diag(Sigma_k);%m*1
                  temp_Sigma_k_ii = sqrt(diag(Sigma_k));
                  phi_ik =   Sigma_k(:,k).*normpdf((b_k(j,:) -mu_k(j,:))'./temp_Sigma_k_ii)./temp_Sigma_k_ii;%m*1
                  for i=1:M
                        % need c.i^(k) and Sigma.i^(k)
                        cik =  b_k(j,:)'  - mu_k(j,:)' -Sigma_k(:,i).* temp_ci ; cik(i)=[]; % m-1*1
                        sigmaik =  (Sigma_k - Sigma_k(:,i)*Sigma_k(i,:)./Sigma_k(i,i)) ;  sigmaik(i,:)=[]; sigmaik(:,i)=[]; % m-1*1
                        Phi_ik =   qsimvnv(200*M,sigmaik ,-inf.*ones(M-1,1),cik ) ; %1*1
                        ETI(j) =  ETI(j) +  phi_ik(i)*Phi_ik;
                  end
            end
        end
    end
    
end