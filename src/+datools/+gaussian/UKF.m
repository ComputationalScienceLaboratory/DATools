classdef UKF < datools.gaussian.GaussianFilter

    properties
        Alpha
        Beta
        Kappa
        Name = "Unscented Kalman Filter";
    end

    methods

        function obj = UKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            addParameter(p, 'Alpha', 1);
            addParameter(p, 'Beta', 2);
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;
            obj.Alpha = s.Alpha;
            obj.Beta = s.Beta;

            kept = p.Unmatched;

            p = inputParser;
            n = s.Model.NumVars;
            addParameter(p, 'InitialState', zeros(n, 1));
            addParameter(p, 'InitialCovariance', eye(n, n));
            addParameter(p, 'Kappa', 3 - n);
            parse(p, kept);

            s = p.Results;

            obj.MeanEstimate = s.InitialState;
            obj.CovarianceEstimate = s.InitialCovariance;
            obj.Kappa = s.Kappa;
        end

        function forecast(obj)
            alpha = obj.Alpha;
            kappa = obj.Kappa;
            beta = obj.Beta;
            n = size(obj.MeanEstimate, 1);
            lambda = (alpha^2)*(n + kappa) - n;

            % set the local mean and covariance
            Xmu = obj.MeanEstimate;
            P =  obj.CovarianceEstimate;
            % sample sigma points
            sqBt = sqrt(n + lambda)*sqrtm(P);
            Xi = [Xmu, Xmu + sqBt, Xmu - sqBt];

            %FXi = zeros(size(Xi));

            %for i = 1:size(Xi, 2)
            %    [~ , yend] = obj.Model.solve([], Xi(:, i));
            %
            %    FXi(:, i) = yend;
            %end

            [~ , FXi] = obj.Model.solve([], Xi);

            Wm = [lambda/(lambda + n), (1/(2*(n + lambda)))*ones(1, 2*n)];
            Wc = [lambda/(lambda + n) + (1 - alpha^2 + beta), ...
                (1/(2*(n + lambda)))*ones(1, 2*n)];

            FXimu = FXi*(Wm.');
            A = (Xi - Xmu);

            P = (Wc.*A)*(A.') + obj.Model.Uncertainty.Covariance;
            obj.MeanEstimate = FXimu;

            % recompute P
            P = (P + P.')/2;
            obj.CovarianceEstimate = P;
        end

        function analysis(obj, obs)

            y = obs.Uncertainty.Mean;
            R = obs.Uncertainty.Covariance;

            alpha = obj.Alpha;
            kappa = obj.Kappa;
            beta = obj.Beta;
            n = size(obj.MeanEstimate, 1);
            lambda = (alpha^2)*(n + kappa) - n;

            % set the local mean and covariance
            Xmu = obj.MeanEstimate;
            P =  obj.CovarianceEstimate;
            % sample sigma points
            sqBt = sqrt(n + lambda)*sqrtm(P);
            Xi = [Xmu, Xmu + sqBt, Xmu - sqBt];
            HXi = obs.observeWithoutError(Xi);

            Wm = [lambda/(lambda + n), (1/(2*(n + lambda)))*ones(1, 2*n)];
            Wc = [lambda/(lambda + n) + (1 - alpha^2 + beta), ...
                (1/(2*(n + lambda)))*ones(1, 2*n)];

            HXimu = HXi*(Wm.');
            A = (Xi - Xmu);
            Z = (HXi - HXimu);

            S = ((Wc.*Z)*Z.' + R);
            S = (S + S.')/2;

            d = (HXimu - y);
            PHt = (Wc.*A)*(Z.');
            obj.MeanEstimate = Xmu - PHt*(S\d);
            % recompute P
            P = (P - PHt*(S\(PHt.')));
            P = (P + P.')/2;
            obj.CovarianceEstimate = P;
        end
    end


end