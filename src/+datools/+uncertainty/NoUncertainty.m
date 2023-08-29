classdef NoUncertainty < handle

    properties
        Mean = 0
        Covariance = 0
    end

    methods

        function x = sample(~, N)
            x = zeros(1, N);
        end

        function xp = addError(~, x)
            xp = x;
        end

    end

end