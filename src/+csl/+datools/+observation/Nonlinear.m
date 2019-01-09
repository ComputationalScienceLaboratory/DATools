classdef Nonlinear < csl.datools.observation.Observation
    
    properties
        F
    end
    
    methods
        
        function obj = Nonlinear(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'F', @(~, x) x); %, @(f) nargin(f) == 2
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj@csl.datools.observation.Observation(nvars, p.Unmatched);
            
            obj.F = s.F;
            
        end
        
        function y = observeWithoutError(obj, t, x)
            y = obj.F(t, x);
        end
        
    end
    
end