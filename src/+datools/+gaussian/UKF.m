classdef UKF < datools.gaussian.Gaussian

    properties
        Alpha
        Beta
        Kappa
    end

    methods

        function obj = UKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            %addParameter(p, 'ModelUncertainty', datools.uncertainty.NoUncertainty);
            addParameter(p, 'Alpha', 1);
            addParameter(p, 'Beta', 2);
            parse(p, varargin{:});

            s = p.Results;

            modelUncertainty = s.Model.Uncertainty;

            obj.Model = s.Model;
            obj.ModelError = modelUncertainty;
            obj.Alpha = s.Alpha;
            obj.Beta = s.Beta;

            kept = p.Unmatched;

            p = inputParser;
            n = s.Model.NumVars;
            addParameter(p, 'Observation', datools.observation.Observation(n));
            addParameter(p, 'InitialState', zeros(n, 1));
            addParameter(p, 'InitialCovariance', eye(n, n));
            addParameter(p, 'Kappa', 3 - n);
            parse(p, kept);

            s = p.Results;

            obj.Observation = s.Observation;
            obj.State = s.InitialState;
            obj.Covariance = s.InitialCovariance;
            obj.Kappa = s.Kappa;
        end

        function forecast(obj)
            alpha = obj.Alpha;
            kappa = obj.Kappa;
            beta = obj.Beta;
            n = size(obj.State, 1);
            lambda = (alpha^2)*(n + kappa) - n;

            % set the local mean and covariance
            Xmu = obj.State;
            P =  obj.Covariance;
            % sample sigma points
            sqBt = sqrt(n + lambda)*sqrtm(P);
            Xi = [Xmu, Xmu + sqBt, Xmu - sqBt];

            FXi = zeros(size(Xi));

            for i = 1:size(Xi, 2)
                [~ , yend] = obj.Model.solve([], Xi(:, i));

                FXi(:, i) = yend;
            end

            Wm = [lambda/(lambda + n), (1/(2*(n + lambda)))*ones(1, 2*n)];
            Wc = [lambda/(lambda + n) + (1 - alpha^2 + beta), ...
                (1/(2*(n + lambda)))*ones(1, 2*n)];

            FXimu = FXi*(Wm.');
            A = (Xi - Xmu);

            P = (Wc.*A)*(A.') + obj.ModelError.Covariance;
            obj.State = FXimu;

            % recompute P
            P = (P + P.')/2;
            obj.Covariance = P;
        end

        function analysis(obj, obs)

            y = obs.Mean;
            R = obs.Covariance;

            alpha = obj.Alpha;
            kappa = obj.Kappa;
            beta = obj.Beta;
            n = size(obj.State, 1);
            lambda = (alpha^2)*(n + kappa) - n;

            % set the local mean and covariance
            Xmu = obj.State;
            P =  obj.Covariance;
            % sample sigma points
            sqBt = sqrt(n + lambda)*sqrtm(P);
            Xi = [Xmu, Xmu + sqBt, Xmu - sqBt];
            HXi = obj.Observation.observeWithoutError(Xi);

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
            obj.State = Xmu - PHt*(S\d);
            % recompute P
            P = (P - PHt*(S\(PHt.')));
            P = (P + P.')/2;
            obj.Covariance = P;
        end
    end


end