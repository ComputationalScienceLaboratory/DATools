classdef ThreeDVar < handle
    %

    properties
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        ModelError % type err
        Observation % type of obervation
        State % current state values for all the states
        B % background covariance
        OptAlg % lbfgs or newton
    end

    properties (Dependent)
        BestEstimate % Current estimate of the particles/ensembles
    end

    properties (Access=protected)
        BDecomposition
    end

    methods
        function obj = ThreeDVar(varargin)
            %ENF   The constructor initializes the properties/attributes
            %
            %   OBJ = ENF(VARARGIN) accepts variable length argument list
            %   VARARGIN and updates the properties/attributes of the object
            %   (OBJ) of this class or a derived class

            optAlgValFcn = @(x) (strcmp(x, 'lbfgs') || strcmp(x, 'newton'));

            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            addParameter(p, 'ModelError', datools.error.Error);
            addParameter(p, 'InitialState', []);
            addParameter(p, 'BackgroundCovariance', []);
            addParameter(p, 'OptimizationType', 'lbfgs', optAlgValFcn);
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;
            obj.ModelError = s.ModelError;
            obj.State = s.InitialState;
            obj.B = s.BackgroundCovariance;
            obj.BDecomposition = decomposition(obj.B, 'chol');
            obj.OptAlg = s.OptimizationType;

            kept = p.Unmatched;

            p = inputParser;
            addParameter(p, 'Observation', datools.observation.Observation(s.Model.NumVars));
            parse(p, kept);

            s = p.Results;

            obj.Observation = s.Observation;
        end

        function forecast(obj)
            %FORECAST Method to propagate the model forward in time
            %
            %   FORECAST(OBJ) propoagates the model one step in time
            %   using a user defined time integration method

            [~ , yend] = obj.Model.solve([], obj.State);

            obj.State = obj.ModelError.adderr(obj.Model.TimeSpan(end), yend);


        end

        function analysis(obj, R, y)

            % A constrained problem like double pendulum will need a
            % constraint coming in from the object.

            dB = obj.BDecomposition;

            if strcmp(class(R), "decomposition")
                dR = R;
            else
                dR = decomposition(R, 'chol');
            end

            xb = obj.State;
            obs = obj.Observation;

            t = obj.Model.TimeSpan(end);

            H = @(x) obs.observeWithoutError(t, x);
            Hadjoint = @(x) obs.linearization(t, x).';

            J = @(x) cost(x, xb, dB, y, dR, H, Hadjoint);

            if strcmp(obj.OptAlg, 'lbfgs')
                opts = optimoptions('fmincon','Display','none', ...
                    'SpecifyObjectiveGradient',true, ...
                    'HessianApproximation', {'lbfgs', 30});
            else

                hessv = @(Hax, v) dB\v + Hax*(dR\(Hax.'*v));
                opts = optimoptions('fmincon','Display','none', ...
                    'Algorithm','trust-region-reflective', ...
                    'SpecifyObjectiveGradient',true, ...
                    'HessianMultiplyFcn', hessv, ...
                    'SubproblemAlgorithm', 'cg');

            end

            xa = fmincon(J, xb, [], [], [], [], [], [], [], opts);

            obj.State = xa;

            function [c, g, hessinfo] = cost(x, xb, dB, y, dR, H, Hadjoint)

                dx = x - xb;
                dy = H(x) - y;
                Hax = Hadjoint(x);

                Bix = dB\dx;
                Riy = dR\dy;

                c = 0.5*((dx.'*Bix) + (dy.'*Riy));
                g = Bix + Hax*Riy;
                hessinfo = Hax;

            end


        end

        function x = get.BestEstimate(obj)
            %GET.BESTESTIMATE Method to estimate of ensemble values
            %
            %   X = GET.BESTESTIMATE(OBJ) uses an in-built getter method,
            %   derived from handle class, to return the best estimate of
            %   the information from the current ensembles of states and
            %   its corresponding weights

            x = obj.State;

        end

    end

end
