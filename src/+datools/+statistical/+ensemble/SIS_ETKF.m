classdef SIS_ETKF < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            
            inflation = obj.Inflation;
            tau = obj.Rejuvenation;
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            w = obj.Weights;
            
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
            
            modelsize = numel(xam);
            
            
            if modelsize > ensN
                xa = xam + Aa*randn(ensN, ensN);
            else
                Pa = Aa*Aa.';
                xa = xam + sqrtm(Pa)*randn(modelsize, ensN);
            end
            
            % SIS
            Hxa = obj.Observation.observeWithoutError(tc, xa);
            t0 = Hxa - y;
            lik = exp(-0.5*sum(t0.*(R\t0), 1));
            
            t0 = xa - xf;
            prop = exp(-0.5*sum(t0.*(Pa\t0), 1));
            
            w = w.*lik.*prop;
            
            % If model error is present, need to calculate the
            % probabilities of evolution. otherwise, it is assumed 1.  
            
            w = w./sum(w);

            P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
            xa = xa + xf*P;
            
            obj.Ensemble = xa;
            obj.Weights = w;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
