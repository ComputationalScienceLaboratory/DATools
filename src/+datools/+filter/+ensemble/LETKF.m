classdef LETKF < datools.filter.ensemble.EnF
    % Localized Ensemble Transform Kalman Filter
    % citation/reference
    properties
        Name = "Localized Ensemble Transform Kalman Filter"
    end

    methods

        function analysis(obj, obs)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)

            inflation = obj.Inflation;

            tc = obj.Model.ODEModel.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;

            y = obs.Uncertainty.Mean;
            R = obs.Uncertainty.Covariance;

            % R = obs.ErrorModel.Covariance;

            % y = obs.Y;

            xfm = mean(xf, 2);
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation * Af / sqrt(ensN-1);
            xf = repmat(xfm, 1, ensN) + Af;

            Hxf = obs.observeWithoutError(xf);
            Hxfm = mean(Hxf, 2);
            HAf = Hxf - repmat(Hxfm, 1, ensN);

            Hi = obs.Indices;

            Aa = zeros(size(Af));
            xam = zeros(size(xfm));

            invR = spdiags(1./diag(R), 0, size(R, 1), size(R, 2));

            for k = 1:numel(xfm)

                if isempty(obj.Localization)
                    C = speye(size(invR));
                else
                    C = obj.Localization(xfm, Hi, k);
                end

                T = sqrtm(speye(ensN)+HAf.'*(C * invR)*HAf);
                Aa(k, :) = Af(k, :) / T;
                xam(k) = xfm(k) + Aa(k, :) * (HAf / T).' * (C * invR) * (y - Hxfm);

            end

            xa = sqrt(ensN-1) * Aa + repmat(xam, 1, ensN);
            obj.Ensemble = xa;
            obj.Model.update(0, obj.MeanEstimate);

        end

    end

end
