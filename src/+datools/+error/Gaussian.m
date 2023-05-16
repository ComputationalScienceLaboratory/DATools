classdef Gaussian < datools.error.Uncertainty
    %GAUSSIAN Gaussian child class for adding gaussina noise with a user
    %   defined deviation and bias

    properties
        Covariance % Covariance for the noise
        Bias % Bias
    end

    properties (Dependent)
        CovarianceSqrt % Covariance square root
    end

    methods

        function obj = Gaussian(varargin)
            p = inputParser;
            addOptional(p, 'Covariance', 1);
            addOptional(p, 'Bias', 0);
            parse(p, varargin{:});

            s = p.Results;

            obj.Covariance = s.Covariance;
            obj.Bias = s.Bias;
        end


        function error = sample(obj, ~, x)
            error = obj.Bias + obj.CovarianceSqrt * randn(size(x), 'like', x);
        end

        
        function updateCovariance(obj, newCovariance)
            obj.Covariance = newCovariance;
        end

        function sqrt = get.CovarianceSqrt(obj)
            % find other methods/decomposition like 'cholesky', etc.
            sqrt = sqrtm(obj.Covariance);
        end

    end

end
