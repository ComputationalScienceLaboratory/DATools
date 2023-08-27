classdef (Abstract) DABase < handle


    properties (Abstract)
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        Observation % type of obervation
        MeanEstimate
        CovarianceEstimate
        Name
    end


    methods (Abstract)
        forecast(obj)
        analysis(obj, obs)
    end

end
