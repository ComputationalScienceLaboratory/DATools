classdef Gaussian < datools.error.Error
    
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
            A = obj.CovarianceSqrt*randn(numel(x), 1);
            %Am = mean(A, 2);
            %A = A - repmat(Am, 1, size(A, 2));
            xp = x + A + obj.Bias;
        end
        
    end
    
end
