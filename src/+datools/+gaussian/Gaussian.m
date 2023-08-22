classdef Gaussian < datools.DAmethod

    properties
        Model 
        ModelError 
        Observation 
        State
        Covariance
    end

    properties (Dependent)
        BestEstimate % Current estimate of the particles/ensembles
    end

    methods
        function forecast(obj)

            [~ , yend] = obj.Model.solve([], obj.State);

            obj.State = obj.ModelError.adderr(obj.Model.TimeSpan(end), yend);

        end

        function x = get.BestEstimate(obj)

            x = obj.State;

        end
    end
end
