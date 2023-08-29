classdef EnGMF < datools.statistical.ensemble.EnF

    methods

        function analysis(obj, R, y)

            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            n = size(xf, 1);

            beta = ((4/(n + 2))^(2/(n + 4)))*((ensN)^(-2/(n + 4)));

            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = sqrt(beta) * Af/sqrt(ensN - 1);

            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            HAf = sqrt(beta)*HAf/sqrt(ensN - 1);

            dS = decomposition((HAf * HAf.') + R, 'chol');

            xtilde = xf - (Af*HAf.')*(dS\(Hxf - y));
            Aa = xtilde - mean(xtilde, 2);
            Aa = sqrt(beta)*Aa/sqrt(ensN - 1);

            t0 = (Hxf - y);

            as = (-0.5 * sum(t0.*(dS \ t0), 1)).';
            m = max(as);
            w = exp(as-(m + log(sum(exp(as-m))))).';
            w = w/sum(w);

            defensivefactor = 1;
            w = defensivefactor*w + (1 - defensivefactor)/ensN;

            % resample
            xa = xf;
            for  k = 1:ensN
                ind = find(rand <= cumsum(w), 1, 'first');
                xa(:, k) = xtilde(:, ind) + Aa*randn(ensN, 1);
            end

            obj.Ensemble = xa;

            obj.Weights = ones(ensN, 1) / ensN;

        end

    end

end
