classdef LETKF < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            %% Group sections
            inflation = obj.Inflation;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation*Af/sqrt(ensN - 1);
            xf = repmat(xfm, 1, ensN) + Af;
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            
            Hi = obj.Observation.Indicies;
            
            Aa = zeros(size(Af));
            xam = zeros(size(xfm));
            
            invR = spdiags(1./diag(R), 0, size(R, 1), size(R, 2));
            
            for k = 1:numel(xfm)
                
                if isempty(obj.Localization)
                    C = speye(size(invR));
                else
                    C = obj.Localization(tc, xfm, Hi, k);
                end
                
                T = sqrtm(speye(ensN) + HAf.'*(C*invR)*HAf);
                Aa(k, :) = Af(k, :)/T;
                xam(k) = xfm(k) + Aa(k, :)*(HAf/T).'*(C*invR)*(y - Hxfm);
                
            end
            
            xa = sqrt(ensN - 1)*Aa + repmat(xam, 1, ensN);
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
