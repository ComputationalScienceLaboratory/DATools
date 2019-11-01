classdef Model < handle
    
    properties
        ODEModel
        Solver
        SynthError
    end
    
    properties (Dependent)
        State
        NumVars
        TimeSpan
        DistanceFunction
    end
    
    methods
        function obj = Model(varargin)
            
            
            p = inputParser;
            addRequired(p, 'ODEModel');%, @(x) isa(x, 'csl.odetestproblems.Problem'));
            addRequired(p, 'Solver', @(x) nargin(x) == 3);
            addParameter(p, 'SynthError', csl.datools.error.Error);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.ODEModel = s.ODEModel;
            obj.Solver = s.Solver;
            obj.SynthError = s.SynthError;
            
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
                [t, y] = obj.Solver(@(t, y) obj.ODEModel.Rhs.F(t, y, params), tspan, y0);
            else
                [t, y] = obj.Solver(obj.ODEModel.Rhs.F, tspan, y0);
            end
            
            time = t(end) - t(1);
            yend = y(end, :).';
            
        end
        
        function update(obj, time, y0)
            
           obj.ODEModel.Y0 = obj.SynthError.adderr(obj.ODEModel.TimeSpan(end) + time, y0);
           obj.ODEModel.TimeSpan = obj.ODEModel.TimeSpan + time;
           
        end
        
        function state = get.State(obj)
            
            state = obj.ODEModel.Y0;
            
        end
        
        function numvars = get.NumVars(obj)
            
            numvars = obj.ODEModel.NumVars;
            
        end
        
        function timespan = get.TimeSpan(obj)
            
            timespan = obj.ODEModel.TimeSpan;
            
        end
        
        function distfn = get.DistanceFunction(obj)
            
            distfn = obj.ODEModel.DistanceFunction;
            
        end
    end
    
end

