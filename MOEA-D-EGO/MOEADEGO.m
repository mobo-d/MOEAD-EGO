classdef MOEADEGO < ALGORITHM
% <multi> <real> <expensive>
% MOEA/D with efficient global optimization
% q    ---   5 --- Batch size (# of function evaluations at each generation)
% maxIter --- 50 --- The maximum number of iterations in inner optimization
 


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

% This function is written by Liang Zhao
% https://github.com/mobo-d/MOEAD-EGO

    methods
        function main(Algorithm,Problem)
           %% Parameter setting
            [q,maxIter] = Algorithm.ParameterSet(5,50);
            C0 = 11*Problem.D-1; % number of initial design points
            
           %% Initial hyperparameters for GP
            GPModels = cell(1,Problem.M);   theta = cell(1,Problem.M);
            theta(:) = {(C0 ^ (-1./C0)).*ones(1,Problem.D)};
      
           %% Generate initial design using OLHS or other DOE methods
            Decs = srgtsESEAOLHSdesign(C0,Problem.D,50,5);% OLHS (npoints, ndv,maxiter, maxstalliter);
            % Decs = UniformPoint(C0,Problem.D,'Latin'); %  LSH in PlatEMO
            D_pop = Problem.Evaluation(repmat(Problem.upper-Problem.lower,C0,1).*Decs+repmat(Problem.lower,C0,1));
            [FrontNo,~] = NDSort(D_pop.objs,1); % find non-dominated solutions
          
            %% Optimization
            while Algorithm.NotTerminated(D_pop(FrontNo==1))
                D_decs = D_pop.decs; D_objs = D_pop.objs;
              %% Step 1: Bulid GP model for each objective function  
                for i = 1 : Problem.M
                    GPModels{i}= Dacefit(D_decs,D_objs(:,i),'regpoly0','corrgauss',theta{i},1e-6*ones(1,Problem.D),20*ones(1,Problem.D));
                    theta{i}   = GPModels{i}.theta;
                end 
              
              %% Step 2: Maximize ETI using MOEA/D and select q candidate points
                SelectDecs     = Opt_ETI(Problem,GPModels,D_decs,D_objs,maxIter,q);
  
              %% Step 3： Aggregate data
                D_pop = [D_pop,Problem.Evaluation(SelectDecs)];
                [FrontNo,~] = NDSort(D_pop.objs,1);
            end
        end
    end
end