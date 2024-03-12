classdef (Abstract) Uncertainty < handle

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
