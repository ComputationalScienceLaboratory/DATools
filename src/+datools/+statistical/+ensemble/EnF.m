classdef EnF < handle

    properties
        Model % ODE
        ModelError % type err
        Observation
        Ensemble
        Weights
        Inflation
        Rejuvenation
        Localization
        Parallel
        RankHistogram
        RankValue
        ResamplingThreshold
    end

    properties (Dependent)
        BestEstimate
        NumEnsemble
    end

    methods (Abstract)

        analysis(obj, R, y)

    end


    methods
        function obj = EnF(varargin)

            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            addParameter(p, 'ModelError', datools.error.Error);
            addParameter(p, 'NumEnsemble', 1);
            addParameter(p, 'Inflation', 1);
            addParameter(p, 'Rejuvenation', []);
            addParameter(p, 'Localization', []);
            addParameter(p, 'Parallel', false);
            addParameter(p, 'RIPIterations', 0);
            addParameter(p, 'RankHistogram', [])
            addParameter(p, 'ResamplingThreshold', 0.5);
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;
            obj.ModelError = s.ModelError;
            obj.Inflation = s.Inflation;
            obj.Rejuvenation = s.Rejuvenation;
            obj.Localization = s.Localization;
            obj.Parallel = s.Parallel;
            obj.RankHistogram = s.RankHistogram;
            obj.ResamplingThreshold = s.ResamplingThreshold;
            ensN = s.NumEnsemble;
            if ~isempty(obj.RankHistogram)
                RankValue = zeros(length(obj.RankHistogram), ensN+2);
            end

            kept = p.Unmatched;

            p = inputParser;
            addParameter(p, 'Observation', datools.observation.Observation(s.Model.NumVars));
            addParameter(p, 'EnsembleGenerator', @(x) randn(s.Model.NumVars, x));
            parse(p, kept);

            s = p.Results;

            obj.Ensemble = s.EnsembleGenerator(ensN);
            obj.Observation = s.Observation;
            obj.Weights = ones(ensN, 1) / ensN;

        end

        function forecast(obj)

            times = zeros(obj.NumEnsemble, 1);

            if obj.Parallel
                rhs = obj.Model.ODEModel.Rhs.F;
                tspan = obj.Model.TimeSpan;
                solver = obj.Model.Solver;
                ens = obj.Ensemble;
                ensN = obj.NumEnsemble;

                parfor ensi = 1:ensN

                    [t, y] = solver(rhs, tspan, ens(:, ensi));

                    time = t(end) - t(1);
                    yend = y(end, :).';

                    ens(:, ensi) = yend;
                    times(ensi) = time;

                end

                for ensi = 1:ensN
                    obj.Ensemble(:, ensi) = obj.ModelError.adderr(obj.Model.TimeSpan(end), ens(:, ensi));
                end

            else
                for ensi = 1:obj.NumEnsemble
                    [time, yend] = obj.Model.solve([], obj.Ensemble(:, ensi));

                    obj.Ensemble(:, ensi) = obj.ModelError.adderr(obj.Model.TimeSpan(end), yend);
                    times(ensi) = time;
                end
            end

            obj.Model.update(mean(times), obj.BestEstimate);

        end


        function ensN = get.NumEnsemble(obj)

            ensN = size(obj.Ensemble, 2);

        end


        function x = get.BestEstimate(obj)

            x = obj.Ensemble * obj.Weights;

        end


        function setMean(obj, xam)

            X = obj.Ensemble;
            ensN = size(X, 2);
            xm = mean(X, 2);
            A = (X - repmat(xm, 1, ensN));
            X = repmat(xam, 1, ensN) + A;
            obj.Ensemble = X;

        end

        function scaleAnomalies(obj, scale)

            X = obj.Ensemble;
            ensN = size(X, 2);
            xm = mean(X, 2);
            A = scale * (X - repmat(xm, 1, ensN));
            X = repmat(xm, 1, ensN) + A;
            obj.Ensemble = X;

        end

        function rejuvenate(obj, tau, Xf)

            X = obj.Ensemble;
            [n, ensN] = size(X);

            if n > ensN + 2
                A = (Xf - mean(Xf, 2)) / sqrt(ensN-1);
                vs = sqrt(sum(A.^2, 2));

                Xi = sqrt(tau) * vs .* randn(n, ensN);
                Xi = Xi - mean(Xi, 2);
    
                X = X + Xi;
            else
                P = sqrt(tau/(ensN - 1)) * (eye(ensN) - ones(ensN) / ensN) * randn(ensN) * (eye(ensN) - ones(ensN) / ensN);
                X = X + Xf * P;
            end

            obj.Ensemble = X;

        end

    end

end
