classdef ETKF < csl.dataassimilation.statistical.ensemble.EnKF
    
    methods
        function obj = ETKF(problem,solver,H,ensembleSize,E, Qest)
            obj = obj@csl.dataassimilation.statistical.ensemble.EnKF(problem,...
                solver,H,ensembleSize,E, Qest);
        end
        function analysis(obj, observation, R)
           
            
            E = obj.CurrentBestGuess;
            H = obj.H;
            ensN = obj.EnsembleSize;
            y = observation;
            
            % THIS CAN BE CHANGED
            U = eye(ensN);
            
            oneM = ones(1, ensN);
            sem1 = sqrt(ensN - 1);

            xbar = mean(E, 2);
            
            
            X = (E - (xbar * oneM))/sem1;
            
            Z = E(H, :);
            
            ybar = mean(Z, 2);
            
            %sqrtR = sqrtm(R);
            isqrtR = sqrtm(inv(R));
            
            %S = (sqrtR\(Z - (ybar * oneM)))/sem1;
            %delta = sqrtR\(y - ybar);
            
            S = (isqrtR*(Z - (ybar * oneM)))/sem1;
            delta = isqrtR*(y - ybar);
            
            Ti = (eye(ensN) + S'*S);
            T = inv(Ti);
            
            w = Ti\(S')*delta;
            
            %w = T * (S')*delta;
            
            %E = xbar*oneM + X*(w*oneM+sem1*sqrtm(Ti)\U);
            
            E = xbar*oneM + X*(w*oneM+sem1*sqrtm(T)*U);
       
            obj.CurrentBestGuess = E;
            
            % The mean of the analysis is our current understanding of the state
            % of the system.
            obj.Problem.Y0 = mean(E,2);
            
            
        end
    end
    
end