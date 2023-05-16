classdef DEnKF < datools.statistical.ensemble.EnF

    methods

        function analysis(obj)

            inflation = obj.Inflation;

            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            R = obj.Observation.Covariance;

            xfm = mean(xf, 2);

            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation * Af;

            xf = repmat(xfm, 1, ensN) + Af;

            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);

            HAf = Hxf - repmat(Hxfm, 1, ensN);

            % tapering

            if isempty(obj.Localization)
                rhoHt = ones(size(Af, 1), size(HAf, 1));
                HrhoHt = ones(size(HAf, 1), size(HAf, 1));
            else
                H = obj.Observation.linearization(tc, xfm);
                rhoHt = obj.Localization(tc, xfm, H);
                HrhoHt = H * rhoHt;
            end

            PfHt = rhoHt .* ((1 / (ensN - 1)) * (Af * (HAf.')));
            HPfHt = HrhoHt .* ((1 / (ensN - 1)) * (HAf * (HAf.')));

            S = HPfHt + R;

            dS = decomposition(S, 'chol');

            Aa = Af - 0.5 * (PfHt * (dS \ HAf));

            d = obj.Observation.Y - Hxfm;
            xam = xfm + (PfHt * (dS \ d));

            obj.Ensemble = repmat(xam, 1, ensN) + Aa;
            obj.Model.update(0, obj.BestEstimate);

        end

    end

end
