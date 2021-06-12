classdef SIS_EnKF_res < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)
            
            inflation = obj.Inflation;
            tau = obj.Rejuvenation;
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            w = obj.Weights;
            
            % EnKF
            xfm = mean(xf, 2);
            Af = inflation/sqrt(ensN - 1)*(xf - xfm);
            xf = xfm + Af*sqrt(ensN - 1);
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = 1/sqrt(ensN - 1)*(Hxf - Hxfm);
            
            if isempty(obj.Localization)
                rhoHt  = ones(size(Af, 1), size(HAf, 1));
                HrhoHt = ones(size(HAf, 1), size(HAf, 1));
            else
                H = obj.Observation.linearization(tc, xfm);
                rhoHt = obj.Localization(tc, xfm, H);
                HrhoHt = H*rhoHt;
            end
            
            PfHt = rhoHt.*((1/(ensN - 1))*(Af*(HAf.')));
            HPfHt = HrhoHt.*((1/(ensN - 1))*(HAf*(HAf.')));
            
            S = HPfHt + R;
            dS = decomposition(S, 'chol');
            
            inn = y - Hxfm;
            
            u = xf + PfHt*(dS\inn);
            xa = u + PfHt*(dS\(sqrtm(R)*randn(size(y, 1), ensN)));
            
            t0 = xa - u;
            V = (1/(ensN - 1))*(t0*t0.');
            % To deal with low rank V
            V = V + trace(V)/3*(eye(3));
            prop = exp(-0.5*sum(t0.*(V\t0), 1));
            
            Hxa = obj.Observation.observeWithoutError(tc, xa);
            t0 = Hxa - y;
            lik = exp(-0.5*sum(t0.*(R\t0), 1));
            
            % If model error is present, need to calculate the
            % probabilities of evolution. otherwise, it is assumed 1.
            
            w = w.*lik./prop;
            w = w./sum(w);
            
            if 1/sum(w.^2) < ensN
                what = cumsum(w);
                a = rand(1, ensN);
                for i = 1:ensN
                    ind = find(a(i) < what, 1);
                    xa(:, i) = xf(:, ind);
                end
                w = (1/ensN)*ones(1, ensN);
            end
            
            P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
            xa = xa + xf*P;
            
            obj.Ensemble = xa;
            obj.Weights = w;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
