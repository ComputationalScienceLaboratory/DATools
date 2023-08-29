classdef ETKF < datools.statistical.ensemble.EnF

    methods

        function analysis(obj, R, y)

            inflation = obj.Inflation;

            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;

            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation*Af/sqrt(ensN - 1);
            xf = repmat(xfm, 1, ensN) + Af*sqrt(ensN - 1);

            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            HAf = HAf/sqrt(ensN - 1);

            dS = decomposition((HAf*HAf.') + R, 'chol');
            dR = decomposition(R, 'chol');

            T = sqrtm(eye(ensN) - (HAf.'*(dS\HAf)));

            Aa = Af*T;
            xam = xfm + ((Aa*(HAf*T).')*(dR\(y - Hxfm)));
            xa = repmat(xam, 1, ensN) + sqrt(ensN - 1)*Aa;

            obj.Ensemble = xa;

            obj.Weights = ones(ensN, 1)/ensN;

        end

    end

end
