classdef DEnKF < csl.datools.DAmethod
    
    properties
        EnsErr
        EnsembleSize
        LocalisationFunction
        InflationValue
        PreviousAnalysis
        PreviousTimeSpan
        CurrentObsErr
        CurrentSyntheticObs
        CurrentObs
        Qestimate
        IsQError
        Iteration
        Parallel
        RipIts
    end
    
    methods
        function obj = DEnKF(problem, ... % An ODE test problem 
                solver,              ... % An ODE solver
                H,                   ... % A vector type of linear observation operator
                ensembleSize,        ... % The size of the ensemble to use
                E,                   ... % How to perturb the ensemble
                Qest,                ... % Estimate of model error
                locF,                ... % Optional argument for the localisation function
                infV, parallel, ripits)                    % Inflation Value
            
            obj = obj@csl.datools.DAmethod(problem, solver, H);
            
            obj.EnsembleSize = ensembleSize;
            obj.EnsErr = E;
            obj.Qestimate = Qest;
            
            if nargin < 9
                parallel = true;
            end
            
            obj.Parallel = parallel;
            
            ensN  = obj.EnsembleSize;
            numVars = obj.Problem.NumVars;
            
            y0 = obj.Problem.Y0;
            
            obj.Iteration = 0;
            
            % We build the ensemble. This does require the use of Q. Q can be an
            % approximation.
            if isa(E, 'function_handle')
                Xf = E(ensN);
            else
                Xf = mvnrnd(y0, E, ensN).';
            end
            
            obj.CurrentBestGuess = Xf;
            
            if nargin < 7
                % We will set our localisation function to be a do nothing function
                obj.LocalisationFunction = @(~) ones(numVars);
            else
                obj.LocalisationFunction = locF;
            end
            
            if nargin < 8
                obj.InflationValue = 1;
            else
                obj.InflationValue = infV;
            end
            
            if nargin < 9
                obj.RipIts = [1, 1];
            else
                obj.RipIts = ripits;
            end
            
            if normest(Qest)
                obj.IsQError = true;
            else
                obj.IsQError = false;
            end
            
        end
        function forecast(obj)
            obj.PreviousAnalysis = obj.CurrentBestGuess;
            obj.PreviousTimeSpan = obj.Problem.TimeSpan;
            
            ensN = obj.EnsembleSize;
            Xas = obj.CurrentBestGuess;
            solver = obj.Solver;
            f = obj.Problem.F;
            tspan = obj.Problem.TimeSpan;
            isqerr = obj.IsQError;
            Q = obj.Qestimate;
            if obj.Parallel
                parfor i = 1:ensN
                    Xa =  Xas(:, i);
                    [~, Xf] = solver(f, tspan, Xa);
                    Xf = Xf.';
                    if isqerr
                        Xas(:, i) =  mvnrnd(Xf(:, end), Q).';
                    else
                        Xas(:,i) =  Xf(:,end);
                    end
                    
                    % obj.Problem.DrawFrame(obj.Problem.TimeSpan(1), Xf(:,end), []);
                end
            else
                for i = 1:ensN
                    Xa =  Xas(:, i);
                    [~, Xf] = solver(f, tspan, Xa);
                    Xf = Xf.';
                    if isqerr
                        Xas(:, i) =  mvnrnd(Xf(:, end), Q).';
                    else
                        Xas(:,i) =  Xf(:,end);
                    end
                    
                    % obj.Problem.DrawFrame(obj.Problem.TimeSpan(1), Xf(:,end), []);
                end
            end
            obj.CurrentBestGuess = Xas;
            
            obj.Problem.Y0 = mean(obj.CurrentBestGuess, 2);
            tend = obj.Problem.TimeSpan(end);
            deltat = diff(obj.Problem.TimeSpan);
            obj.Problem.TimeSpan = [tend, tend + deltat];
        end
        function analysis(obj, observation, R)
            obj.CurrentObsErr = R;
            obj.CurrentObs = observation;
           
            H = obj.H;
            ensN = obj.EnsembleSize;
            y = observation;
            
            if obj.Iteration < obj.RipIts(2) 
                ripits = obj.RipIts(1);
            else
                ripits = 1;
            end
            
            for ri = 1:ripits
                Xf =  obj.CurrentBestGuess;
                Xfm = mean(Xf, 2);
                
                if ri == 1
                    Af = Xf - Xfm*ones(1, ensN);
                    Af = obj.InflationValue * Af;
                    Xf = Xfm*ones(1, ensN) + Af;
                    obj.CurrentBestGuess = Xf;
                end
                
                Af = Xf - Xfm*ones(1, ensN);
                
                PfHt = obj.getForecastCovarience(Af);
                
                % S matrix from Law et al.
                S = PfHt(H, :) + R;
                
                % Kalman gain matrix
                K = PfHt/S;
                
                % Hey, analysis step of EnKF
                d = (y - Xfm(H, :));
                Xam = Xfm + K*d;
                
                Aa = Af - 1/2*K*Af(H, :);
                
                %Aa = obj.InflationValue * Aa;
                
                Xa = Xam*ones(1, ensN) + Aa;
                
                obj.CurrentBestGuess = Xa;
                
                
                
                % The mean of the analysis is our current understanding of the state
                % of the system.
                obj.Problem.Y0 = Xam;
                
                
                %obj.Problem.DrawFrame(obj.Problem.TimeSpan(1), obj.Problem.Y0, []);
            end
            
            obj.Iteration = obj.Iteration + 1;
        end
        
        function PfHt = getForecastCovarience(obj, Af)
            ensN = obj.EnsembleSize;

            H = obj.H;
            
            PfHt = 1/(ensN - 1) * (Af*(Af(H, :)'));
            
            % Localise
            rhoHt = obj.getLocalisationMatrix(H);
            
            % Schur product
            PfHt = rhoHt .* PfHt;
            
            %Pf = obj.InflationValue * Pf;
        end
        function rho = getLocalisationMatrix(obj, H)
            rho = obj.LocalisationFunction(obj, H);
        end
    end
    
end