classdef Indexed < datools.observation.Observation
    
    properties
        Indicies
    end
    
    methods
        
        function obj = Indexed(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Indicies', 1);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj@datools.observation.Observation(nvars, p.Unmatched);
            
            obj.Indicies = s.Indicies;
            
        end
        
        function y = observeWithoutError(obj, ~, x)
            y = x(obj.Indicies, :);
        end
        
        function H = linearization(obj, ~, ~)
            I = speye(obj.NumVars);
            H = I(obj.Indicies, :);
        end
        
    end
    
end