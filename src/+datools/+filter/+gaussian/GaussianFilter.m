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

    methods
        function forecast(obj)

            [~ , yend] = obj.Model.solve([], obj.State);

            obj.MeanEstimate = obj.Model.Uncertainty.adderr(obj.Model.TimeSpan(end), yend);

        end

    end
end
