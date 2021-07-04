classdef SIR < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            
            tau = obj.Rejuvenation;
            
            tc = obj.Model.TimeSpan(1);
            
            dR = decomposition(R, 'chol');
            
            xf = obj.Ensemble;
            xa = xf;
            
            ensN = obj.NumEnsemble;            
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            t0 = Hxf - y;
              
            % more efficient way of calculating weights
            as = (-0.5*sum(t0.*(dR\t0), 1)).';
            m = max(as);
            w = exp(as - (m + log(sum(exp(as - m)))));
            
            what = cumsum(w);
            a = rand(1, ensN);
            
            for i = 1:ensN
                ind = find(a(i) < what, 1);
                xa(:, i) = xf(:, ind);
            end
            
            n = size(xa, 1);
            
            if n < ensN + 2
                Af = (xf - mean(xf))/sqrt(ensN -1);
                
                Xi = sqrt(tau)*randn(n, ensN); Xi = Xi - mean(Xi, 2);
                
                vs = sqrt(sum(Af.^2, 2));
                
                Xi = sqrt(tau)*vs.*rand(n, ensN);
                Xi = Xi - mean(Xi, 2);

                xa = xa + Xi;
            else
                P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
                xa = xa + xf*P;
            end
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
            obj.Weights = ones(ensN, 1)/ensN;
            
        end
        
    end
    
end
