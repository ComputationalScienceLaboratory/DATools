classdef ETPF < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)

            tau = obj.Rejuvenation;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            AeqT = [kron(speye(ensN), ones(1, ensN)); kron(ones(ensN, 1), speye(ensN)).'];
            lbT = zeros((ensN)*ensN, 1);
            optsT = optimoptions('linprog', 'Display', 'off');
           
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            
            xdist = zeros(ensN, ensN);
            w = zeros(ensN, 1);
            
            for i = 1:ensN
                xtemp = xf - repmat(xf(:, i), 1, ensN);
                xdist(i, :) = vecnorm(xtemp).^2;
                
                inn = y - Hxf(:, i);
                w(i) = exp(-0.5*inn'*(R\inn));
            end
            
            w = w/sum(w);
            
            beqT = [ones(ensN, 1)/ensN; w];
            f = xdist(:);
            Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
            Tx = ensN*reshape(Tx, ensN, ensN);
            
            P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
            xa = xf*(Tx + P);
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end