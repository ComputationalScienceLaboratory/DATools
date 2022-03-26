classdef EnF < handle
    % This is the base class for all statistical methods
    % Derive from this class and overload methods/functions as required
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
        % this is the method that will be overloaded by child classes
        analysis(obj, R, y)

    end


    methods
        function obj = EnF(varargin)
            % EnF  The constructor initializes the properties/attributes
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
            % forecast Method to propagate the model forward in time and
            %
            % INPUT PARAMETERS:
            % obj : reference to the object of this class/ derived class
            %
            % OUTPUT PARAMETERS:
            % NULL

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
            % get.NumEnsemble   Method to get the number of ensembles (a getter method).
            %                   A getter method is used for dependent properties
            %
            % INPUT PARAMETERS:
            % obj : reference to the object of this class/derived class
            %
            % OUTPUT PARAMETERS:
            % ensN : number of ensemble members

            ensN = size(obj.Ensemble, 2);

        end


        function x = get.BestEstimate(obj)
            % Method to get the best estimates of Ensemble values. A getter
            % method is used to derive the value of a dependent property.
            % Current estimate is made by multiplying weight os each
            % ensembles with its current value.
            %
            % INPUT PARAMETERS:
            % obj : reference to the object of this class/derived class
            %
            % OUTPUT PARAMETERS:
            % x : current ensemble values

            x = obj.Ensemble * obj.Weights;

        end


        function setMean(obj, xam)
            % setMean   Method to set the mean of the ensembles, if required
            %
            % INPUT PARAMETERS:
            % obj : reference to the object of this class/derived class
            %
            % OUTPUT PARAMETERS:
            % Null

            X = obj.Ensemble;
            ensN = size(X, 2);
            xm = mean(X, 2);
            A = (X - repmat(xm, 1, ensN));
            X = repmat(xam, 1, ensN) + A;
            obj.Ensemble = X;

        end

        function scaleAnomalies(obj, scale)
            % scaleAnomalies  Method to scale the anomalies of the ensembles, if required
            %
            % INPUT PARAMETERS:
            % obj : reference to the object of this class/derived class
            % scale : scaling factor
            %
            % OUTPUT PARAMETERS:
            % Null

            X = obj.Ensemble;
            ensN = size(X, 2);
            xm = mean(X, 2);
            A = scale * (X - repmat(xm, 1, ensN));
            X = repmat(xm, 1, ensN) + A;
            obj.Ensemble = X;

        end

        function rejuvenate(obj, tau, Xf)
            % rejuvenate  To reduce particle degeneracy, rejuvenation is equivalent to
            %             adding random noise (in the form of random combination of
            %             background anomalies) to the transformation matrix. Give
            %             appropriate reference to notes and ETPF (to be done).
            %
            % INPUT PARAMETERS:
            % obj : reference to the object of this class/derived class
            % tau : rejuvenation bandwidth
            % xf : forecast/background information
            %
            % OUTPUT PARAMETERS:
            % Null

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
