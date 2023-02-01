classdef ETKF < datools.statistical.ensemble.EnF

    methods

        function analysis(obj)

            inflation = obj.Inflation;

            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            R = obj.Observation.R;

            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation * Af/sqrt(ensN - 1);
            xf = repmat(xfm, 1, ensN) + Af*sqrt(ensN - 1);

            Hxf = obj.ObservationOperator.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            HAf = HAf/sqrt(ensN - 1);

            dS = decomposition((HAf * HAf.')+ R, 'chol');
            dR = decomposition(R, 'chol');

            temp = dS \ HAf;
            T = sqrtm(eye(ensN)-(HAf.' * temp));

            Aa = Af * T;
            xam = xfm + ((Aa * (HAf * T).') * (dR \ (obj.Observation.Y - Hxfm)));
            xa = sqrt(ensN-1) .* Aa + repmat(xam, 1, ensN);

            obj.Ensemble = xa;

            obj.Weights = ones(ensN, 1) / ensN;

        end

    end

end
