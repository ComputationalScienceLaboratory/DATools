classdef LETKF < csl.dataassimilation.DAmethod
    properties
        EnsErr
        EnsembleSize
        LocalisationFunction
        PreviousAnalysis
        PreviousTimeSpan
        CurrentObsErr
        CurrentObs
        Qestimate
    end
    
    methods
        function obj = LETKF(problem, ... % An ODE test problem 
                solver,              ... % An ODE solver
                H,                   ... % A vector type of linear observation operator
                ensembleSize,        ... % The size of the ensemble to use
                E,                   ... % How to perturb the ensemble
                Qest,                ... % Estimate of model error
                locF)
            
            obj = obj@csl.dataassimilation.DAmethod(problem, solver, H);
            
            obj.EnsembleSize = ensembleSize;
            obj.EnsErr = E;
            obj.Qestimate = Qest;
            
            ensN  = obj.EnsembleSize;
            numVars = obj.Problem.NumVars;
            
            y0 = obj.Problem.Y0;
            
            % We build the ensemble. This does require the use of Q. Q can be an
            % approximation.
            Xf = mvnrnd(y0, E, ensN)';
            
            obj.CurrentBestGuess = Xf;
            
            if nargin < 7
                % We will set our localisation function to be a do nothing function
                obj.LocalisationFunction = @(~) ones(numVars);
            else
                obj.LocalisationFunction = locF;
            end
        end
        function forecast(obj)
            obj.PreviousAnalysis = obj.CurrentBestGuess;
            obj.PreviousTimeSpan = obj.Problem.TimeSpan;
            
            ensN = obj.EnsembleSize;
            for i = 1:ensN
                Xa =  obj.CurrentBestGuess(:, i);
                [t, Xf] = obj.Solver(obj.Problem.F, obj.Problem.TimeSpan, Xa);
                Xf = Xf';
                obj.CurrentBestGuess(:, i) =  mvnrnd(Xf(:, end), obj.Qestimate)';
                %obj.CurrentBestGuess(:,i) =  Xf(:,end);
            end
            
            obj.Problem.Y0 = mean(obj.CurrentBestGuess, 2);
            deltat = diff(obj.Problem.TimeSpan);
            obj.Problem.TimeSpan = [t(end), t(end) + deltat];
        end
        function analysis(obj, observation, R)
            obj.CurrentObsErr = R;
            obj.CurrentObs = observation;
            
            ensN = obj.EnsembleSize;
            
            y = observation;
            %obsN = numel(y);
            
            H = obj.H;
            
            % our current best guess is the forecast that we made before.
            Xf =  obj.CurrentBestGuess;
            
            Xfm = mean(Xf, 2);
            
            % Reich and Cotter
            Af = Xf - Xfm*ones(1, ensN);
            
            HAf = Af(H, :);
            
            Xa = zeros(size(Xf));

            for k = 1:numel(Xfm)
                %
                
                Ctilde = obj.generateCTilde(k);
                
                HCtildeH = Ctilde(H, H);
                
                HCtildeHRi = (HCtildeH/R);
                
                Sk2i = eye(ensN, ensN) + 1/(ensN - 1) * (HAf.') * (HCtildeHRi) * HAf;
                
                Ski = sqrtm(Sk2i);
                
                Aak = Af(k, :)/Ski;
                
                wk = 1/ensN - 1/(ensN - 1) * (Sk2i\(HAf.')) * HCtildeHRi * (Xfm(H, :) - y);
              
                
                Xakm = Xf(k, :)*wk;
                
                
                Xa(k, :) = Xakm +Aak;
            end
            
            obj.CurrentBestGuess = Xa;
            
            % The mean of the analysis is our current understanding of the state
            % of the system.
            obj.Problem.Y0 = mean(Xa, 2);
        end

        function Ctilde = generateCTilde(obj, k)
            Ctilde = obj.LocalisationFunction(obj, k);
        end
    end
    
end

