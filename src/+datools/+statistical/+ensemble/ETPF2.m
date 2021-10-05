classdef ETPF2 < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)

            tau = obj.Rejuvenation;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            dR = decomposition(R, 'chol');

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
            
            
            % Delta
            W = diag(w);
            Bric = Tx - w*ones(1, ensN);
            Aric = ensN*(W - w*w.') - Bric*Bric.';
            
            rs = @(D) reshape(D, ensN, ensN);
            urs = @(D) reshape(D, ensN*ensN, 1);
            f = @(t, D) urs(-Bric * rs(D) - rs(D)*(Bric.') + Aric - rs(D)*rs(D));
            D = zeros(ensN, ensN);
            D = rs(datools.utils.dpETPF2(f, D(:)));
            
            xa = xf*(Tx + D);
            
            obj.Ensemble = xa;
            obj.Weights = ones(ensN, 1)/ensN;
            
            obj.rejuvenate(tau, xf);
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
