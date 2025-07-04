classdef FETKF < datools.filter.ensemble.EnF
    % Fancy Ensemble Transform Kalman Filter
    % citation/reference
    properties
        B
        Bsqrt
        BsqrtD
        SurrogateEnsN
        Laplace
        Gamma
        Name = 'Fancy Ensemble Transform Kalman Filter' % Please check
        Type = "Ensemble"
    end

    methods
        function obj = FETKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'B', 1);
            addParameter(p, 'Bsqrt', 1)
            addParameter(p, 'SurrogateEnsembleSize', 2);
            addParameter(p, 'Laplace', false);

            parse(p, varargin{8:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.filter.ensemble.EnF(varargin{1:7}, kept);

            obj.B = s.B;
            obj.Bsqrt = s.Bsqrt;
            obj.BsqrtD = decomposition(obj.Bsqrt, 'chol');

            obj.SurrogateEnsN = s.SurrogateEnsembleSize;
            obj.Laplace = s.Laplace;

        end

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
            
            R = obs.Uncertainty.Covariance;

            xfm = mean(xf, 2);
            n = numel(xfm);


            Af = xf - repmat(xfm, 1, ensN);

            Af = inflation * Af;

            Bs = obj.Bsqrt;
            BsD = obj.BsqrtD;

            s = svd(BsD \ Af/sqrt(ensN-1));
            trBPf = sum(s.^2);
            trBPf2 = sum(s.^4);

            ensNN = ensN - 1;

            if isempty(obj.Gamma)
                % RBLW
                gamma = min((((ensNN - 2) / ensNN) * trBPf2 + trBPf^2)/((ensNN + 2) * (trBPf2 - (1 / n) * trBPf^2)), 1);
                % OAS
                %gamma = max(min(( (1 - 2/n)*trBPf2 + trBPf^2 )/( (ensNN + 1 - 2/n)*(trBPf2 - (trBPf^2)/n) ), 1), 0);
            else
                gamma = obj.Gamma;
            end
            mu = trBPf / n;

            ensNW = obj.SurrogateEnsN;

            xf = repmat(xfm, 1, ensN) + Af;

            Hxf = obs.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);

            HAf = Hxf - repmat(Hxfm, 1, ensN);

            Afs = sqrt(1-gamma) * Af / sqrt(ensN-1);

            w = sqrt(mu) * Bs * randn(n, ensNW);
            w = w - repmat(mean(w, 2), 1, ensNW);
            ws = sqrt(gamma) * w / sqrt(ensNW-1);

            Hw = obs.observeWithoutError(tc, w+repmat(xfm, 1, ensNW));
            Hwm = mean(Hw, 2);
            HAw = Hw - repmat(Hwm, 1, ensNW);

            ZA = sqrt(1-gamma) * HAf / sqrt(ensN-1);
            Zw = sqrt(gamma) * HAw / sqrt(ensNW-1);

            S = ZA * ZA' + Zw * Zw' + R;

            ZZ = [ZA, Zw];

            TT = (eye(ensN+ensNW) - ZZ.' / S * ZZ);

            d = obs.Y - Hxfm;

            AAa = [Afs, ws] * real(sqrtm(TT));
            Aa = sqrt(ensN-1) * AAa(:, 1:ensN) / sqrt(1-gamma);

            xam = xfm + ([Afs, ws] * TT * ZZ.') / R * d;

            obj.Ensemble = repmat(xam, 1, ensN) + Aa;

            obj.Model.update(0, obj.BestEstimate);

        end

    end

end
