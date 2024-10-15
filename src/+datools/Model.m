classdef Model < handle
    %MODEL this is the base model class
    %   Handles all the data pertaining to a given model
    %   User can use OTP model or define their own
    %   User may define their choice of time integration method (or MATLODE)
    %   Model errors can be incorporated using the SynthError

    properties
        ODEModel % refactor as this may not be an ode or otp problem
        Solver % refactor because a test problem can come with its own solver
        SolverTLM % the model operator is linearized in this case
        Uncertainty % the model noise
    end

    properties (Dependent)
        State % current states of the Model
        NumVars % number of state variables
        DistanceFunction % distance function for localization
        TimeStamp % timestamp for this curent state
    end

    methods
        function obj = Model(varargin)

            p = inputParser;
            addRequired(p, 'ODEModel'); %, @(x) isa(x, 'csl.odetestproblems.Problem'));
            addRequired(p, 'Solver', @(x) nargin(x) == 3);
            addParameter(p, 'Uncertainty', datools.uncertainty.NoUncertainty);
            addParameter(p, 'SolverTLM', @(f, t, y, J, lambda) datools.utils.rk4tlm(f, t, y, J, lambda, 10));
            parse(p, varargin{:});

            s = p.Results;

            obj.ODEModel = s.ODEModel;
            obj.Solver = s.Solver;
            obj.Uncertainty = s.Uncertainty;
            obj.SolverTLM = s.SolverTLM;

        end

        function evolve(obj, time)

            tstart = obj.ODEModel.TimeSpan(1);

            if nargin < 2
                time = diff(obj.ODEModel.TimeSpan);
            end

            tend = tstart + time;

            tspan = [tstart, tend];

            [time, yend] = obj.solve(tspan);

            obj.update(time, yend);

        end

        function [time, yend] = solve(obj, tspan, y0, params)

            if nargin < 2 || isempty(tspan)
                tspan = obj.ODEModel.TimeSpan;
            end

            if nargin < 3 || isempty(y0)
                y0 = obj.ODEModel.Y0;
            end

            if nargin < 4 || isempty(params)
                params = [];
            end

            if ~isempty(params)
                [t, yend] = obj.Solver(@(t, y) obj.ODEModel.RHS.F(t, y, params), tspan, y0);
            else
                % [t, yend] = obj.Solver(obj.ODEModel.RHS.F, tspan, y0);
                [t, yend] = obj.Solver(obj.ODEModel.F, tspan, y0);
            end

            time = t(end) - t(1);
            % yend = yend(end, :).';

        end

        function [time, yend, lambda] = solveWithTLM(obj, tspan, y0, lambda, params)

            if nargin < 2 || isempty(tspan)
                tspan = obj.ODEModel.TimeSpan;
            end

            if nargin < 3 || isempty(y0)
                y0 = obj.ODEModel.Y0;
            end

            if nargin < 4 || isempty(lambda)
                lambda = eye(size(y0, 1));
            end

            if nargin < 5 || isempty(params)
                params = [];
            end

            if ~isempty(params)
                % [t, yend] = obj.Solver(@(t, y) obj.ODEModel.RHS.F(t, y, params), tspan, y0);
                [t, yend] = obj.Solver(@(t, y) obj.ODEModel.F(t, y, params), tspan, y0);
            else
                J = obj.ODEModel.RHS.Jacobian;
                if isempty(J)
                    % J = @(t, y) otp.utils.derivatives.jacobian( ...
                    %     obj.ODEModel.RHS.F, t, y, 'FD');
                    J = @(t, y) otp.utils.derivatives.jacobian( ...
                        obj.ODEModel.F, t, y, 'FD');
                end

                % [t, yend, lambda] = obj.SolverTLM(obj.ODEModel.RHS.F, tspan, y0, J, lambda);
                [t, yend, lambda] = obj.SolverTLM(obj.ODEModel.F, tspan, y0, J, lambda);
            end

            time = t(end) - t(1);
            % yend = yend(end, :).';

        end

        function update(obj, time, y0)

            obj.ODEModel.Y0 = obj.Uncertainty.addError(y0);
            obj.ODEModel.TimeSpan = obj.ODEModel.TimeSpan + time;

        end

        function state = get.State(obj)

            state = obj.ODEModel.Y0;

        end

        function numvars = get.NumVars(obj)

            numvars = obj.ODEModel.NumVars;

        end


        function distfn = get.DistanceFunction(obj)

            distfn = obj.ODEModel.DistanceFunction;

        end

        function t = get.TimeStamp(obj)
            t = obj.ODEModel.TimeSpan(1);
        end
    end

end
