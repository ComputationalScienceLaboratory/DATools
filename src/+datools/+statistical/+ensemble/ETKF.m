classdef ETKF < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            
            inflation = obj.Inflation;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
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
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end