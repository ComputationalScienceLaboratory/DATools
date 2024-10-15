classdef IEKF < datools.filter.gaussian.GaussianFilter
    % Iterated Extended Kalman Filter
    % DOI: https://doi.org/10.1109/9.250476
    properties
        Name = "Iterated Extended Kalman Filter";
        Tolerance
    end

    methods

        function obj = IEKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            addParameter(p, 'Tolerance', 1e-3);
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;
            obj.Tolerance = s.Tolerance;

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

            h = @(x) obs.observeWithoutError(x);
            V = @(x, xb, P) 0.5*(h(x) - y).'*(R\(h(x) - y)) ...
                + 0.5*(x - xb).'*(P\(x - xb));

            Xc = X;
            Vcur = V(Xc, X, P);

            gnorm = inf;
            linesearchdone = false;
            while gnorm > obj.Tolerance
                H = obs.linearization(Xc);

                S = (H*P*H.' + R);
                S = (S + S.')/2;

                deltax = X - Xc + P*(H.'*(S\(y - h(Xc) ...
                    - H*(X - Xc))));
                gnorm = norm(deltax);

                % perform naive line search
                alpha = 1;

                while true
                    Xcnew = Xc + alpha*deltax;
                    Vnew = V(Xcnew, X, P);

                    if alpha < eps
                        linesearchdone = true;
                        break
                    end

                    if Vnew < Vcur
                        Xc = Xcnew;
                        Vcur = Vnew;
                        break;
                    else
                        alpha = alpha/(1.5);
                    end
                end

                if linesearchdone
                    break;
                end
            end

            % update the covariance
            H = obs.linearization(Xc);

            S = (H*P*H.' + R);
            S = (S + S.')/2;

            PHt = P*H.';

            obj.MeanEstimate = Xc;
            % recompute P
            P = (P - PHt*(S\(PHt.')));
            P = (P + P.')/2;
            obj.CovarianceEstimate = P;
        end
    end


end