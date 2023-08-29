classdef EKF < datools.gaussian.GaussianFilter

    properties
        Name = "Extended Kalman Filter";
    end

    methods

        function obj = EKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;

            kept = p.Unmatched;

            p = inputParser;
            n = s.Model.NumVars;
            addParameter(p, 'InitialState', zeros(n, 1));
            addParameter(p, 'InitialCovariance', eye(n, n));
            parse(p, kept);

            s = p.Results;

            obj.MeanEstimate = s.InitialState;
            obj.CovarianceEstimate = s.InitialCovariance;
        end

        function forecast(obj)
            n = size(obj.MeanEstimate, 1);

            % set the local mean and covariance
            X = obj.MeanEstimate;
            P =  obj.CovarianceEstimate;


            [~ , FX, F] = obj.Model.solveWithTLM([], X, eye(n));

            P = F*P*F.' + obj.Model.Uncertainty.Covariance;
            obj.MeanEstimate = FX;

            % recompute P
            P = (P + P.')/2;
            obj.CovarianceEstimate = P;
        end

        function analysis(obj, obs)
            y = obs.Uncertainty.Mean;
            R = obs.Uncertainty.Covariance;

            % set the local mean and covariance
            X = obj.MeanEstimate;
            P =  obj.CovarianceEstimate;

            H = obs.linearization(X);
            
            S = H*P*H.' + R;

            PHt = P*H.';

            HX = obs.observeWithoutError(X);

            d = HX - y;

            obj.MeanEstimate = X - PHt*(S\d);
            % recompute P
            P = (P - PHt*(S\(PHt.')));
            P = (P + P.')/2;
            obj.CovarianceEstimate = P;
        end
    end


end