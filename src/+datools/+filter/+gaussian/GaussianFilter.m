classdef GaussianFilter < datools.DABase
    %GAUSSIANFILTER This is the base class for all gaussian based methods
    %   Derive from this class and implement methods/functions as required
    %   Deriving from handle base class allows an object of this class to be
    %   passed by reference.
    properties
        Model
        MeanEstimate
        CovarianceEstimate
    end

    properties (Abstract)
        Name % name of the filter
        Type % type of filter(Ensemble/Gaussian/Variational)
    end

    methods (Abstract)
        % A method that will be implemented by child  classes to make
        % approximate inference on ensembles of states by combining
        % prior forecast/background data with noisy observations
        analysis(obj, obs)
    end

    methods
        function forecast(obj)

            [~, yend] = obj.Model.solve([], obj.State);

            obj.MeanEstimate = obj.Model.Uncertainty.adderr(obj.Model.TimeSpan(end), yend);

        end

    end
end
