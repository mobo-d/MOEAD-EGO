classdef MOEADEGO < ALGORITHM
% <multi> <real> <expensive>
% MOEA/D-EGO: adaptive adjust z* and direction vectors
% batch_size    ---   5 --- number of true function evaluations per iteration
% n_init   ---   0 --- number of initial samples. A value of 0 indicates that n_init=11d-1. 

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

% This function was written by Liang Zhao (liazhao5-c@my.cityu.edu.hk).
 

    methods
        function main(Algorithm,Problem)
           %% Parameter setting
            [batch_size, n_init] = Algorithm.ParameterSet(5,0); 
            % parameters for MOP
            D = Problem.D; M = Problem.M;
            xlower = Problem.lower;
            xupper = Problem.upper;
            % number of initial samples
            if n_init == 0, n_init = 11*D-1;end 
            % Initial hyperparameters for GP
            theta = repmat({(n_init ^ (-1 ./ n_init)) .* ones(1, D)}, 1, M);
            
            %% Generate initial design using LHS or other DOE methods
            x_lhs = lhsdesign(n_init, D,'criterion','maximin','iterations',1000);
            x_init = xlower +  (xupper - xlower).*x_lhs;  
            Archive = Problem.Evaluation(x_init);     
            % find non-dominated solutions
            [FrontNo,~] = NDSort(Archive.objs,1); 
          
            %% Optimization
            while Algorithm.NotTerminated(Archive(FrontNo==1))
                train_x = Archive.decs; train_y = Archive.objs;
              %% Bulid GP model for each objective function  
                GPModels = cell(1,M);
                for i = 1 : M
                    GPModels{i}= Dacefit(train_x,train_y(:,i),'regpoly0','corrgauss',theta{i},1e-6*ones(1,D),20*ones(1,D));
                    theta{i}   = GPModels{i}.theta;
                end 
              
              %% Maximize ETI using MOEA/D and select q candidate points
                Batch_size = min(Problem.maxFE - Problem.FE,batch_size); % the total budget is  Problem.maxFE 
                new_x = Opt_ETI(M,D,xlower,xupper,GPModels,Batch_size,train_x,train_y);  
                
 
               %% Expensive Evaluation
                Archive = [Archive,Problem.Evaluation(new_x)];
                [FrontNo,~] = NDSort(Archive.objs,1);
            end
        end
    end
end