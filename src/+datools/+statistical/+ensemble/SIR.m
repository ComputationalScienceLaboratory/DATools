classdef SIR < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            
            % abuse
            tau = obj.Inflation;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            xa = 0*xf;
            
            ensN = obj.NumEnsemble;            
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            t0 = Hxf - y;
            lik = exp(-0.5*sum(t0.*(R\t0)));
            w = lik/sum(lik);
            
            what = cumsum(w);
            a = rand(1, ensN);
            
            for i = 1:ensN
                ind = find(a(i) < what, 1);
                xa(:, i) = xf(:, ind);
            end
            
            P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
            xa = xa + xf*P;
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
