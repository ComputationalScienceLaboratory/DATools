classdef EnF < handle
    % This is the base class for all statistical methods
    % Derive from this class and implement methods/functions as required
    % Deriving from handle base class allows an object of this class to be
    % passed by reference.

    properties
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        ModelError % type err
        Observation % type of obervation
        Ensemble % current ensemble values for all the states
        Weights % Weight of each particle/ensemble
        Inflation % inflation constant
        Rejuvenation % rejuvenation constant
        Localization % A boolean if localization needs to be used
        Parallel % A boolean if parallel threads are to be implemented
        RankHistogram % State variables for which RH is needed
        RankValue % store the RH values for each state variables required
        ResamplingThreshold % threshold below which resampling needs to be done
    end

    properties (Dependent)
        BestEstimate % Current estimate of the particles/ensembles
        NumEnsemble % Number of ensemble members
    end

    methods (Abstract)
        % A method that will be implemented by child  classes to make
        % approximate inference on ensembles of states by combining
        % prior forecast/background data with noisy observations
        analysis(obj, R, y)
    end


    methods
        function obj = EnF(varargin)
            %ENF   The constructor initializes the properties/attributes
            %
            %   OBJ = ENF(VARARGIN) accepts variable length argument list
            %   VARARGIN and updates the properties/attributes of the object
            %   (OBJ) of this class or a derived class

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
            %FORECAST Method to propagate the model forward in time
            %
            %   FORECAST(OBJ) propoagates the model one step in time
            %   using a user defined time integration method

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
            %GET.NUMENSEMBLES   Method to get the number of ensembles
            %
            %   ENSN = GET.NUMENSEMBLES(OBJ) uses an in-built getter method,
            %   derived from handle class, to get the current number of
            %   ensembles of states for objects of this class or a derived
            %   class

            ensN = size(obj.Ensemble, 2);

        end


        function x = get.BestEstimate(obj)
            %GET.BESTESTIMATE Method to estimate of ensemble values
            %
            %   X = GET.BESTESTIMATE(OBJ) uses an in-built getter method,
            %   derived from handle class, to return the best estimate of
            %   the information from the current ensembles of states and
            %   its corresponding weights

            x = obj.Ensemble * obj.Weights;

        end


        function setMean(obj, xam)
            %SETMEAN   Method to set the mean of the ensembles, if required
            %
            %   SETMEAN(OBJ, XAM) sets the mean of the ensembles of states
            %   of the object of this class or a derived class to XAM

            X = obj.Ensemble;
            ensN = size(X, 2);
            xm = mean(X, 2);
            A = (X - repmat(xm, 1, ensN));
            X = repmat(xam, 1, ensN) + A;
            obj.Ensemble = X;

        end

        function scaleAnomalies(obj, scale)
            % SCALEANOMALIES  Method to scale the anomalies of the ensembles
            %
            %   SCALEANOMALIES(OBJ, SCALE) scales the unbiased current
            %   ensembles of state using the scaling factor SCALE

            X = obj.Ensemble;
            ensN = size(X, 2);
            xm = mean(X, 2);
            A = scale * (X - repmat(xm, 1, ensN));
            X = repmat(xm, 1, ensN) + A;
            obj.Ensemble = X;

        end

        function rejuvenate(obj, tau, Xf)
            % REJUVENATE  To reduce particle degeneracy
            %
            %   REJUVENATE(OBJ, TAU, XF) adds a random noise in the form of
            %   random combination of background anomalies of the ensembles
            %   of states (using rejuvenation bandwidth TAU) to the
            %   transformation matrix. This matrix is used to rejuvenate
            %   the collasping states XF and prevent particle degeneracy

            X = obj.Ensemble;
            [n, ensN] = size(X);

            if n > ensN + 2
                A = (Xf - mean(Xf, 2)) / sqrt(ensN-1);
                vs = sqrt(sum(A.^2, 2));

                Xi = sqrt(tau) * vs .* randn(n, ensN);
                Xi = Xi - mean(Xi, 2);

                X = X + Xi;
            else
                P = sqrt(tau/(ensN - 1)) * (eye(ensN) - ones(ensN) / ensN) ...
                    * randn(ensN) * (eye(ensN) - ones(ensN) / ensN);
                X = X + Xf * P;
            end

            obj.Ensemble = X;

        end

    end

end
