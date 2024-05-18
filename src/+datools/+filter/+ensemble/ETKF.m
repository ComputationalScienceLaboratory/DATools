classdef ETKF < datools.filter.ensemble.EnF
    properties
        Name = "Ensemble Transform Kalman Filter"
    end

    methods

        function analysis(obj, obs)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)

            y = obs.Uncertainty.Mean;
            R = obs.Uncertainty.Covariance;
            
            inflation = obj.Inflation;

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            % inflate the background covariance
            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation*Af/sqrt(ensN - 1);
            xf = repmat(xfm, 1, ensN) + Af*sqrt(ensN - 1);
            
            % observation space
            Hxf = obs.observeWithoutError(xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            HAf = HAf/sqrt(ensN - 1);

            dS = decomposition((HAf * HAf.')+ R, 'chol');

            dR = decomposition(R, 'chol');

            T = sqrtm(eye(ensN) - (HAf.'*(dS\HAf)));

            Aa = Af * T;
            xam = xfm + ((Aa * (HAf * T).') * (dR \ (y - Hxfm)));
            xa = sqrt(ensN-1) .* Aa + repmat(xam, 1, ensN);

            obj.Ensemble = xa;

            obj.Weights = ones(ensN, 1)/ensN;

        end

    end

end
