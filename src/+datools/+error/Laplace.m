classdef Laplace < datools.error.Error
    
    properties
        CovarianceSqrt
        Bias
    end
    
    methods
        
        function obj = Laplace(varargin)
            p = inputParser;
            addOptional(p, 'CovarianceSqrt', 1);
            addOptional(p, 'Bias', 0);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.CovarianceSqrt = s.CovarianceSqrt;
            obj.Bias = s.Bias;
        end
        
        function xp = adderr(obj, ~, x)
            X = obj.CovarianceSqrt*randn(size(x, 1), size(x, 2));
            Z = exprnd(1, 1, size(x, 2));
            Y = sqrt(Z).*X;
            xp = x + Y + obj.Bias;
        end
        
    end
    
end
