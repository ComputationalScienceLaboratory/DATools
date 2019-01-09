classdef Gaussian < csl.datools.error.Error
    
    properties
        Rsqrt
    end
    
    methods
        
        function obj = Gaussian(varargin)
            p = inputParser;
            addOptional(p, 'CovarianceSqrt', 1);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.Rsqrt = s.CovarianceSqrt;
        end
        
        function xp = adderr(obj, ~, x)
            xp = x + obj.Rsqrt*randn(numel(x), 1);
        end
        
    end
    
end
