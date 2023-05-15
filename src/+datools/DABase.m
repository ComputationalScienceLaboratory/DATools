classdef DABase < handle
%DABASE This is the base class for all data assimilation problems
%   provide other info

    properties
        Model                 % Model and solver information
        ModelError            % Model Error
        AddErrorToModel       % boolean if we consider/ignore model error
    end
    
    methods
        function obj = DABase(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'Model', []);
            addParameter(p, 'ModelError', []);
            addParameter(p, 'AddErrorToModel', false);
            
            parse(p, varargin{:});

            s = p.Results;
            
            obj.Model = s.Model;
            obj.ModelError = s.ModelError;
            obj.AddErrorToModel = s.AddErrorToModel;
        end
    end
    
    methods(Abstract)
        forecast(obj, observation)
    end
        
end