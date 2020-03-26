classdef ETPF < datools.statistical.ensemble.EnF
    
    
    methods
        
        function analysis(obj, R, y)
            
            tau = obj.Inflation^2 - 1;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            
            xdist = zeros(ensN, ensN);
            weights = zeros(ensN, 1);
            
            for i = 1:ensN
                xtemp = xf - repmat(xf(:, i), 1, ensN);
                xtemp = xtemp .* xtemp;
                xdist(i, :) = sum(xtemp);
                
                inn = y - Hxf(:, i);
                weights(i) = exp(-0.5*inn'*(R\inn));
            end
            
            
            weights = weights/sum(weights)
            
            T = optimvar('T', ensN, ensN, 'LowerBound', 0);
            prob = optimproblem('Objective', sum(sum(T.*xdist)), 'ObjectiveSense', 'min');
            opts = optimoptions(prob);
            opts.Display = 'off';
%             opts.Algorithm = 'interior-point';
            prob.Constraints.c1 = sum(T) == ones(1, ensN)/ensN;
            prob.Constraints.c2 = sum(T, 2) == weights;
            prob.Constraints.c3 = sum(sum(T)) == 1;
            problem = prob2struct(prob, 'Options', opts);
            [sol, ~, ~, ~] = linprog(problem);
            
            Tx = reshape(sol, ensN, ensN);
            Tx = ensN*Tx;
            
            xa = xf*Tx;
            
            B = tau/(ensN - 1)*(Af*Af.');
            xa = xa + sqrtm(B)*randn(size(xa));
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
