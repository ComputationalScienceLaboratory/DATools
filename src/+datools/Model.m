classdef Model < handle
    %MODEL this is the base model class
    %   Handles all the data pertaining to a given model
    %   User can use OTP model or define their own
    %   User may define their choice of time integration method (or MATLODE)
    %   Model errors can be incorporated using the SynthError

    properties
        ODEModel % The ODE model from OTP or user defined
        Solver % Time Integration Scheme
        SynthError % Model Error
        AddError % boolean to add/not add noise during forward propogation
        TimeSpan % Time Span for Propogating the Model
    end

    properties (Dependent)
        State % Sates of the Model
        NumVars % Number of State Variables
        DistanceFunction % Distance Function for Localization
    end

    methods
        function obj = Model(varargin)

            p = inputParser;
            addRequired(p, 'ODEModel'); %, @(x) isa(x, 'csl.odetestproblems.Problem'));
            addRequired(p, 'Solver', @(x) nargin(x) == 3);
            addParameter(p, 'SynthError', []);
            addParameter(p, 'AddError', false);
            parse(p, varargin{:});

            s = p.Results;

            obj.ODEModel = s.ODEModel;
            obj.Solver = s.Solver;
            obj.SynthError = s.SynthError;
            obj.AddError = s.AddError;
            obj.TimeSpan = obj.ODEModel.TimeSpan;

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
                [t, y] = obj.Solver(@(t, y) obj.ODEModel.F(t, y, params), tspan, y0);
            else
                [t, y] = obj.Solver(obj.ODEModel.F, tspan, y0);
            end

            time = t(end) - t(1);
            yend = y(end, :).';

        end

        function update(obj, time, y0)
            if obj.AddError
                obj.ODEModel.Y0 = obj.SynthError.addError(obj.ODEModel.TimeSpan(end)+time, y0);
                obj.ODEModel.TimeSpan = obj.ODEModel.TimeSpan + time;
            else
                obj.ODEModel.Y0 = obj.SynthError.addNoError(obj.ODEModel.TimeSpan(end)+time, y0);
                obj.ODEModel.TimeSpan = obj.ODEModel.TimeSpan + time;
            end

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
    end

end
