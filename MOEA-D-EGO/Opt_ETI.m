function SelectDecs = Opt_ETI(Problem,model,D_decs,D_objs,q)
% Maximizing N Subproblems and Selecting Batch of Points 
% Expected Tchebycheff Improvement (ETI)

%------------------------------- Reference --------------------------------
% Q. Zhang, W. Liu, E. Tsang, and B. Virginas, Expensive multiobjective
% optimization by MOEA/D with Gaussian process model, IEEE Transactions on
% Evolutionary Computation, 2010, 14(3): 456-474.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao.
% https://github.com/mobo-d/MOEAD-EGO

   %% Generate the initial weight vectors
    % if Problem.M > 3, it is recommended to use incremental lattice design
    % or two-layered simplex-lattice design to generate the weight vectors.
    [W, Problem.N]  = UniformPoint(Problem.N,Problem.M); % simplex-lattice design 
    
    %% Parameter setting for MOEA/D-DE
    delta=0.9; nr = 2; maxIter = 50;

    %% update Z^*
    % zi is set as the minimum values of the predicted means obtained by the GP model of i-th objective
    [~,~,Z] = MOEAD_GP(Problem,maxIter,W,D_objs,model,delta,nr); 
   
    %% Using MOEA/D-DE to Maximize ETI
    [Pop_ETI,PopDec,~,~ ] = MOEAD_ETI(Problem,maxIter,W,D_objs,model,delta,nr,Z);
   
    %% Select the unsimilar candidate solutions and build candidate pool Q
    Q = []; Q_ETI = [];  temp = D_decs;
    for i = 1 : Problem.N
        if min(pdist2(real(PopDec(i,:)),real(temp))) > 1e-5
            if Pop_ETI(i)>0
                Q = [Q;PopDec(i,:)]; Q_ETI = [Q_ETI;Pop_ETI(i)];
                temp = [temp;PopDec(i,:)];
            end
        end
    end
    %% Kmeans cluster the solutions into q clusters and select the solutions with the maximum ETI in each cluster
    Batch_size = min(Problem.maxFE - Problem.FE,q); % the total budget is Problem.maxFE 
    SelectDecs = K_means_Batch_Select(Q,Batch_size,PopDec,Q_ETI) ;

end
 
function  [PopDec,Pop_u,Z] = MOEAD_GP(Problem,maxIter,W,D_objs,model,delta,nr) 
%% using MOEA/D to solve subproblems
    Z       = min(D_objs,[],1);
  
    %% neighbourhood   
    T       = ceil(Problem.N/10); % size of neighbourhood
    B       = pdist2(W,W);
    [~,B]   = sort(B,2);
    B       = B(:,1:T);
     %% the initial population for MOEA/D
    PopDec = repmat(Problem.upper-Problem.lower,Problem.N,1).*UniformPoint(Problem.N,Problem.D,'Latin')+repmat(Problem.lower,Problem.N,1); 
    %% initial population
    [Pop_u,~] = Evaluate(PopDec,model);
     Z       = min(min(Pop_u,[],1),Z);
    for gen = 1 : maxIter-1
       for i = 1 : Problem.N
           if rand < delta
               P = B(i,randperm(size(B,2)));
           else
               P = randperm(Problem.N);
           end
           %% generate an offspring and calculate its predictive mean and s
           OffDec = OperatorDE(Problem,PopDec(i,:),PopDec(P(1),:),PopDec(P(2),:));          
           [Off_u,~]= Evaluate(OffDec,model);  
           
            Z = min(Z,Off_u);        
            g_old = max(abs(Pop_u(P,:) - repmat(Z,length(P),1)).*W(P,:),[],2);
            g_new = max(repmat(abs(Off_u-Z),length(P),1).*W(P,:),[],2);
             
           % Update the solutions in P
           offindex = P(find(g_old>g_new,nr));
           if ~isempty(offindex)
               PopDec(offindex,:) = repmat(OffDec,length(offindex),1); 
               Pop_u(offindex,:) = repmat(Off_u,length(offindex),1);
           end
       end      
    end
end

function  [Pop_ETI,PopDec,Pop_u,Pop_s] = MOEAD_ETI(Problem,maxIter,W,D_objs,model,delta,nr,Z) 

%% using MOEA/D to solve subproblems
    N = size(W,1); % N: # of subproblems
    gmin  = Calculate_gmin(D_objs,Z,W);
    %% neighbourhood   
    T       = ceil(N/10); % size of neighbourhood
    B       = pdist2(W,W);
    [~,B]   = sort(B,2);
    B       = B(:,1:T);
    %% initial population, Generate random population
    PopDec = repmat(Problem.upper-Problem.lower,N,1).*UniformPoint(N,Problem.D,'Latin')+repmat(Problem.lower,N,1); 
    %% GP prediction
    [Pop_u,Pop_s] = Evaluate(PopDec,model);
     Pop_ETI = ETICal(Pop_u,Pop_s,W,gmin,Z);   
    for gen = 1 : maxIter-1
       for i = 1 : N
           if rand < delta
               P = B(i,randperm(size(B,2)));
           else
               P = randperm(N);
           end
          %% generate an offspring and calculate its predictive mean and s
           OffDec = OperatorDE(Problem,PopDec(i,:),PopDec(P(1),:),PopDec(P(2),:));          
           [Off_u,Off_s]= Evaluate(OffDec,model);  
        
           g_new = ETICal(repmat(Off_u,length(P),1),repmat(Off_s,length(P),1),W(P,:),gmin(P),Z);
            
           offindex =  find(Pop_ETI(P)<g_new,nr) ;
           if ~isempty(offindex)
               PopDec(P(offindex),:) = repmat(OffDec,length(offindex),1); 
               Pop_u(P(offindex),:) = repmat(Off_u,length(offindex),1);
               Pop_s(P(offindex),:) = repmat(Off_s,length(offindex),1);
               Pop_ETI(P(offindex)) = g_new(offindex);
           end
       end      
    end
end 
function SelectDecs = K_means_Batch_Select(Q,q,PopDec,Q_ETI) 
     batch_size = min(q,size(Q,1));% in case Q is smaller than Batch size
    
    if batch_size == 0
        Qb = randperm(size(PopDec,1),q);
        SelectDecs = PopDec(Qb,:);
    else
        cindex  = kmeans(Q,batch_size);
        Qb = [];
        for i = 1 : batch_size
            index = find(cindex == i); 
            [~,best] = max(Q_ETI(index));
            Qb = [Qb,index(best)];
        end
        SelectDecs = Q(Qb,:);
    end
    % when Q is smaller than batch size
    if size(SelectDecs,1) < q
        Qb = randperm(size(PopDec,1), q - size(SelectDecs,1));
        SelectDecs = [SelectDecs;PopDec(Qb,:)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  gmin  = Calculate_gmin(Popobjs,Z,W)
 % calcate the gmin using the true samples
 % g(x|lambda) = max{lambda1(f1-z1),lambda2(f2-z2)}
    [~,m] = size(W); % n is the number of sub-problems
    Objs = Popobjs-Z;
    G = W(:,1)*Objs(:,1)';
    for j = 2:m
        G = max(G,W(:,j)*Objs(:,j)');
    end
    gmin = min(G,[],2);
end 

function  ETI = ETICal(u,sigma,W,Gbest,Z)
%     g(x|w) = max{w(f1-z1),w(f2-z2)}  
% calculate the ETI(x|w) at multiple requests, e.g., n  
% u     : n*M  predictive mean
% sigma : n*M  square root of the predictive variance
% W     : n*M  weight vectors 
% Gbest : n*1  
% Z     : 1*M  reference point   

    [n,M] = size(u);
    
    g_mu = W.*(u - repmat(Z,n,1));% n*M
    g_sig = W.*sigma; % n*M
    Method = 0;
    if Method == 0 % Approximation
        ETI = Approx_ETI(g_mu,g_sig,Gbest);
    elseif Method == 1 % Exact formula
        ETI = Cal_ETI(g_mu,g_sig,Gbest);
    end
     
end

function [u,s] = Evaluate(X,model)
% Predict the objective vector of the candidate solutions 
    N = size(X,1); % number of samples
    M = length(model); % number of objectives
    u = zeros(N,M); % predictive mean
    MSE = zeros(N,M); % predictive SME
    if N == 1 
        for j = 1 : M
            [u(:,j),~,MSE(:,j)] = Predictor(X,model{1,j}); % DACE Kriging toolbox
        end
        MSE(MSE<0) = 0;
    else
        for j = 1 : M
            [u(:,j),MSE(:,j)] = Predictor(X,model{1,j}); % DACE Kriging toolbox
        end
        MSE(MSE<0) = 0;
    end
   s = sqrt(MSE);% square root of the predictive variance
end
