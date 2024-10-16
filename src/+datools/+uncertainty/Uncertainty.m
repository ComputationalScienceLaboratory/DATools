classdef (Abstract) Uncertainty < handle
    %UNCERTAINTY This is the base class for all noises. Currently we
    % support the following kind of uncertainty
    % Uniform, Gaussian, Laplace, Gaussian Mixture, Empirical
    % 
    properties (Abstract)
        Mean
        Covariance
    end

    methods
        sample(obj, N);
        addError(obj, x)
        asGaussian(obj)
        asGMM(obj)
        pdf(obj)
        log(obj)
    end

end
