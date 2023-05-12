classdef DABase < handle
%DABASE This is the base class for all data assimilation problems
%   provide other info

    properties
        Model                 % Model and solver information
        ModelError            % Model Error
        AddErrorToModel       % boolean if we consider/ignore model error
        ObservationOperator   % Observation operator used (linear, non-linear, etc)
        Observation           % Curent snapshots of reality
    end
    
    methods
        function obj = DABase(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'Model', []);
            addParameter(p, 'ModelError', []);
            addParameter(p, 'AddErrorToModel', false);
            addParameter(p, 'ObservationOperator', []);
            addParameter(p, 'Observation', []);
            
            parse(p, varargin{:});

            s = p.Results;
            
            obj.Model = s.Model;
            obj.ModelError = s.ModelError;
            obj.AddErrorToModel = s.AddErrorToModel;
            obj.ObservationOperator = s.ObservationOperator;
            obj.Observation = s.Observation;
        end
    end
    
    methods(Abstract)
        forecast(obj)
    end
        
end