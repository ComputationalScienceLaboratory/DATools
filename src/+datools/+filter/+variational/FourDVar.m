classdef FourDVar < datools.filter.variational.Variational
    properties
        Name = 'FourDVar'
        N % observation window
    end

    methods
        function analysis(obj, obs)
            y = obs.Uncertainty.Mean;
            R = obs.Uncertainty.Covariance;
            dR = decomposition(R, 'chol');

            % TODO
        end
    end
end