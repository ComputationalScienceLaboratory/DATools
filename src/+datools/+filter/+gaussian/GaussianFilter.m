classdef GaussianFilter < datools.DABase

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
