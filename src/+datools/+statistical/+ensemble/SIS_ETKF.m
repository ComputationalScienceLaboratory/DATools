classdef SIS_ETKF < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            
            inflation = obj.Inflation;
            tau = obj.Rejuvenation;
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            % ETKF
            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation*Af./sqrt(ensN - 1);
            xf = repmat(xfm, 1, ensN) + Af;
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            
            temp = ((HAf*HAf.') + R)\HAf;
            T = sqrtm(eye(ensN) - (HAf.'*temp));
            
            Aa = Af*T;
            xam = xfm + ((Aa*(HAf*T).')*(R\(y - Hxfm)));
            xa = sqrt(ensN - 1).*Aa + repmat(xam, 1, ensN);
            
            xb = xa;       
            
            % SIS
            Hxb = obj.Observation.observeWithoutError(tc, xb);
            t0 = Hxb - y;
            w = exp(-0.5*sum(t0.*(R\t0), 1));
            w = w/sum(w);
            
            what = cumsum(w);
            a = rand(1, ensN);
            
            for i = 1:ensN
                ind = find(a(i) < what, 1);
                xa(:, i) = xb(:, ind);
            end
            
            P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
            xa = xa + xf*P;
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
