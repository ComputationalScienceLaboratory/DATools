classdef Error < handle
%ERROR base class for all kind of errors
%   Curently supports Gaussian, Laplace, Logistics, Tent and Triangle.

    methods

        function xp = adderr(~, ~, x)
            xp = x;
        end

    end

end
