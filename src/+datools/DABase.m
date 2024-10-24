classdef (Abstract) DABase < handle
    % DABase This is the base class for all data assimilation algorithms
    %   Users need to inherit this class and define the properties and
    %   methods in child/derived class

    properties (Abstract)
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        MeanEstimate % current estimates of the states
        CovarianceEstimate % spread of the states, if needed
        Name % the name of the filter
    end


    methods (Abstract)
        forecast(obj) % forecasts the model 
        analysis(obj, obs) % assimilation algorithm
    end

end
