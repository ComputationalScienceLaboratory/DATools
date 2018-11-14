classdef EnKF < csl.dataassimilation.DAmethod
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
    end
    
    methods
        function obj = EnKF(problem, ... % An ODE test problem 
                solver,              ... % An ODE solver
                H,                   ... % A vector type of linear observation operator
                ensembleSize,        ... % The size of the ensemble to use
                E,                   ... % How to perturb the ensemble
                Qest,                ... % Estimate of model error
                locF,                ... % Optional argument for the localisation function
                infV)                    % Inflation Value
            
            obj = obj@csl.dataassimilation.DAmethod(problem, solver, H);
            
            obj.EnsembleSize = ensembleSize;
            obj.EnsErr = E;
            obj.Qestimate = Qest;
            
            ensN  = obj.EnsembleSize;
            numVars = obj.Problem.NumVars;
            
            y0 = obj.Problem.Y0;
            
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
            for i = 1:ensN
                Xa =  obj.CurrentBestGuess(:, i);
                [t, Xf] = obj.Solver(obj.Problem.F, obj.Problem.TimeSpan, Xa);
                Xf = Xf.';
                if obj.IsQError
                    obj.CurrentBestGuess(:, i) =  mvnrnd(Xf(:, end), obj.Qestimate).';
                else
                    obj.CurrentBestGuess(:,i) =  Xf(:,end);
                end
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
            
            % Synthetic observations          
            Ys = mvnrnd(y, R, ensN).';
            
            obj.CurrentSyntheticObs = Ys;
            
            % our current best guess is the forecast that we made before.
            Xf =  obj.CurrentBestGuess;
            Xfm = mean(Xf, 2);
            
            Af = Xf - Xfm*ones(1, ensN);
            Af = obj.InflationValue * Af;
            Xf = Xfm*ones(1, ensN) + Af;
            
            obj.CurrentBestGuess = Xf;
            
            PfHt = obj.getForecastCovarience(Af);

            % S matrix from Law et al.
            S = PfHt(H, :) + R;
            
            % Kalman gain matrix
            K = PfHt/S;

            % Hey, analysis step of EnKF
            d = (Ys - Xf(H, :));
            Xa = Xf + K*d;
            
            
            
            % WE INFLATE THE ANALYSIS ENSEMBLE DEVIATIONS      
            %Xam = mean(Xa, 2);
            
            %Aa = Xa - Xam*ones(1, ensN);
            
            %Aa = obj.InflationValue*Aa;
            
            
            %Xa = Xam*ones(1, ensN) + Aa;
            %
            %
            obj.CurrentBestGuess = Xa;

            
            % The mean of the analysis is our current understanding of the state
            % of the system.
            obj.Problem.Y0 = mean(Xa, 2);
            
            
            %obj.Problem.DrawFrame(0, obj.Problem.Y0, []);
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

