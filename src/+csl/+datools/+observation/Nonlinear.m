classdef Nonlinear < csl.datools.observation.Observation
    
    properties
        F
        J
    end
    
    methods
        
        function obj = Nonlinear(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'F', @(~, x) x);
            addParameter(p, 'J', @(~, x) speye(numel(x)));
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj@csl.datools.observation.Observation(nvars, p.Unmatched);
            
            obj.F = s.F;
            obj.J = s.J;
            
        end
        
        function y = observeWithoutError(obj, t, x)
            y = obj.F(t, x);
        end
        
        function H = linearization(obj, t, x)
            H = obj.J(t, x);
        end
        
    end
    
end