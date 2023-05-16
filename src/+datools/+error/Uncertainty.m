classdef Uncertainty < handle
%ERROR base class for all kind of errors
%   Curently supports Gaussian, Laplace, Logistics, Tent and Triangle.

    methods

        function xp = addError(obj, ~, x)
            xp = x + obj.sample([], x);
        end

        function xp = addNoError(~, ~, x)
            xp = x;
        end

    end

    methods(Abstract)
        sample(obj)
    end


end
