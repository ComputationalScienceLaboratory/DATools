classdef Gaussian < csl.datools.error.Error
    
    properties
        CovarianceSqrt
        Bias
    end
    
    methods
        
        function obj = Gaussian(varargin)
            p = inputParser;
            addOptional(p, 'CovarianceSqrt', 1);
            addOptional(p, 'Bias', 0);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.CovarianceSqrt = s.CovarianceSqrt;
            obj.Bias = s.Bias;
        end
        
        function xp = adderr(obj, ~, x)
            xp = x + obj.CovarianceSqrt*randn(numel(x), 1) + obj.Bias;
        end
        
    end
    
end
