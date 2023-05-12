classdef Gaussian < datools.error.Uncertainty
%GAUSSIAN Gaussian child class for adding gaussina noise with a user
%   defined deviation and bias

    properties
        CovarianceSqrt    % Covariance(sqrt) for the noise
        Bias              % Bias
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


        function error = sample(obj, ~, x)
            error = obj.Bias + obj.CovarianceSqrt * randn(size(x), 'like', x);
        end

        % function xp = addError(obj, ~, x)
        %     xp = x + obj.sample([], x);
        % end

    end

end
