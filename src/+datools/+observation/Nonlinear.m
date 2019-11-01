classdef Nonlinear < csl.datools.observation.Observation
    
    properties
        F
        %J
        Indicies
    end
    
    methods
        
        function obj = Nonlinear(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'F', @(~, x) x);
            %addParameter(p, 'J', @(~, x) speye(numel(x)));
            addParameter(p, 'Indicies', 1);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj@csl.datools.observation.Observation(nvars, p.Unmatched);
            
            obj.F = s.F;
            %obj.J = s.J;
            obj.Indicies = s.Indicies;
            
        end
        
        function y = observeWithoutError(obj, t, x)
            yfull = obj.F(t, x);
            y = yfull(obj.Indicies, :);
        end
        
        function H = linearization(obj, t, x)
            %H = obj.J(t, x);
            H = 1;
        end
        
    end
    
end