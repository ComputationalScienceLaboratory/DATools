classdef ETPF < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)

            tau = obj.Rejuvenation;
            
            tc = obj.Model.TimeSpan(1);
            
            dR = decomposition(R, 'chol');
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            AeqT = [kron(speye(ensN), ones(1, ensN)); kron(ones(ensN, 1), speye(ensN)).'];
            lbT = zeros((ensN)*ensN, 1);
            optsT = optimoptions('linprog', 'Display', 'off');
           
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            
            xdist = zeros(ensN, ensN);
            
            for i = 1:ensN
                xtemp = xf - repmat(xf(:, i), 1, ensN);
                xdist(i, :) = vecnorm(xtemp).^2;
            end
            
            t0 = Hxf - y;
            
            % more efficient way of calculating weights
            as = (-0.5*sum(t0.*(dR\t0), 1)).';
            m = max(as);
            w = exp(as - (m + log(sum(exp(as - m)))));
            
            beqT = [ones(ensN, 1)/ensN; w];
            f = xdist(:);
            Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
            Tx = ensN*reshape(Tx, ensN, ensN);
            
            xa = xf*Tx;
            
            obj.Ensemble = xa;
            obj.Weights = ones(ensN, 1)/ensN;
            obj.rejuvenate(tau);
            
            obj.Model.update(0, obj.BestEstimate);
            
            
            
        end
        
    end
    
end
