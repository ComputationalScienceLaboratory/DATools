classdef (Abstract) DAmethod < handle


    properties (Abstract)
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        ModelError % type err
        Observation % type of obervation
        BestEstimate
    end


    methods (Abstract)
        forecast(obj)
        % A method that will be implemented by child  classes to make
        % approximate inference on ensembles of states by combining
        % prior forecast/background data with noisy observations
        analysis(obj, R, y)
    end

end
