function ETI = Approx_ETI(g_mu,g_sig,Gbest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the ETI(x|w) at multiple requests, e.g., n  
% g_mu  : n*M  predictive mean
% g_sig : n*M  square root of the predictive variance
% Gbest : n*1   
%------------------------------- Reference --------------------------------
% Q. Zhang, W. Liu, E. Tsang, and B. Virginas, Expensive multiobjective
% optimization by MOEA/D with Gaussian process model, IEEE Transactions on
% Evolutionary Computation, 2010, 14(3): 456-474.
%--------------------------------------------------------------------------
% This function is written by Liang Zhao
% https://github.com/mobo-d/MOEAD-EGO

    [~,M] = size(g_mu);
    g_sig(g_sig<0) = 0; g_sig2 = g_sig.^2; % n*M
     
     % Eq. 18 & Eq. 19 in MOEA/D-EGO
	[y,s2] = App_Max_of_2_Gaussian(g_mu(:,1:2),g_sig2(:,1:2)); % f1 & f2
    
    if M >= 3
        for i = 3 : M
            mu_temp = [y,g_mu(:,i)]; sig2_temp = [s2,g_sig2(:,i)];
            [y,s2] = App_Max_of_2_Gaussian(mu_temp,sig2_temp);
        end
    end
    s = (sqrt(s2));
    temp = ((Gbest-y)./s); % n*1
    ETI = (Gbest-y).*normcdf(temp) + s.*normpdf(temp);
   
end
function [y,s2] = App_Max_of_2_Gaussian(mu,sig2)
% Calculate  Eq. 18 & Eq. 19 in MOEA/D-EGO 
% n requests
% mu is n*2
% sig2 is n*2
    tao = sqrt(sum(sig2,2));  % n*1
    alpha = (mu(:,1)-mu(:,2))./tao;  % n*1
    % Eq. 16 / Eq. 18
    y = mu(:,1).*normcdf(alpha) + mu(:,2).*normcdf(-alpha) + tao.*normpdf(alpha);  % n*1
    % There is a typo in Eq. 17.  See Appendix B of MOEA/D-EGO.
    % It should be $$ +(\mu_1+\mu_2) \tau \varphi(\alpha)$$
    s2 = (mu(:,1).^2 + sig2(:,1)).*normcdf(alpha) + ...
        (mu(:,2).^2 + sig2(:,2)).*normcdf(-alpha) + (sum(mu,2)).*tao.*normpdf(alpha); 
%     s2 = (mu(:,1).^2 + sig2(:,1)).*normcdf(alpha) + ...
%         (mu(:,2).^2 + sig2(:,2)).*normcdf(-alpha) + (sum(mu,2)).*normpdf(alpha); 
    s2 = s2 - y.^2;
    s2(s2<0) = 0;
end