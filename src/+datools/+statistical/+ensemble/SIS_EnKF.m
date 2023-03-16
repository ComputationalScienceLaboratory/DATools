classdef SIS_EnKF < datools.statistical.ensemble.EnF

    methods

        function analysis(obj)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)
            
            inflation = obj.Inflation;
            tau = obj.Rejuvenation;
            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            w = obj.Weights;

            modelsize = size(xf, 1);
            obsize = size(obj.Observation.Y, 1);
            
            R = obj.Observation.Covariance;

            % EnKF
            xfm = mean(xf, 2);
            Af = inflation / sqrt(ensN-1) * (xf - xfm);
            xf = xfm + Af * sqrt(ensN-1);

            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = 1 / sqrt(ensN-1) * (Hxf - Hxfm);

            if isempty(obj.Localization)
                rhoHt = ones(modelsize, obsize);
                HrhoHt = ones(obsize, obsize);
            else
                H = obj.Observation.linearization(tc, xfm);
                rhoHt = obj.Localization(tc, xfm, H);
                HrhoHt = H * rhoHt;
            end

            PfHt = rhoHt .* ((1 / (ensN - 1)) * (Af * (HAf.')));
            HPfHt = HrhoHt .* ((1 / (ensN - 1)) * (HAf * (HAf.')));

            S = HPfHt + R;
            dS = decomposition(S, 'chol');

            u = xf + PfHt * (dS \ (obj.Observation.Y - Hxfm));
            t0 = PfHt * (dS \ (sqrtm(R) * randn(obsize, ensN)));
            xa = t0 + u;

            V = (1 / (ensN - 1)) * (t0 * t0.');

            % To deal with low rank V when obsize < modelsize
            V = V + trace(V) / modelsize * (eye(modelsize));

            prop = exp(-0.5*sum(t0.*(V \ t0), 1)).';

            Hxa = obj.Observation.observeWithoutError(tc, xa);
            t0 = Hxa - obj.Observation.Y;
            lik = exp(-0.5*sum(t0.*(R \ t0), 1)).';

            % If model error is present, need to calculate the probabilities of evolution
            % otherwise, it is assumed 1.


            w = w .* lik ./ prop;
            w = w ./ sum(w);

            %             1/sum(w.^2)
            % Resample is number of particles is low. Resampling at each
            % step is better, but we can choose a threshold to resample.
            if 1 / sum(w.^2) < ensN / 2
                what = cumsum(w);
                a = rand(1, ensN);
                for i = 1:ensN
                    ind = find(a(i) < what, 1);
                    xa(:, i) = xf(:, ind);
                end
                w = (1 / ensN) * ones(ensN, 1);
            end


            obj.Ensemble = xa;
            obj.Weights = w;
            obj.rejuvenate(tau, xf);

            obj.Model.update(0, obj.BestEstimate);

        end

    end

end
