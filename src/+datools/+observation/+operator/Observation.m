classdef Observation < handle
%OBSERVATION Base class for every observation operator
%   other details 

    properties
        Y             % Current Observation of States (maybe sparse)
        Covariance    % Observation Error Covariance
        NumVars      % number of state variables
        ErrorModel   % observation error model
        NumObs        % Number of Observations(Make sure > 0)
    end

    methods

        function obj = Observation(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Y', []);
            addRequired(p, 'NumVars');
            addParameter(p, 'ErrorModel', []);
            validationFunction = @(x) (x>0);
            addParameter(p, 'NumObs', 1, validationFunction);

            parse(p, varargin{:});

            s = p.Results;

            obj.NumVars = s.NumVars;
            obj.ErrorModel = s.ErrorModel;
            obj.Y = s.Y;
            obj.NumObs = s.NumObs;

            unMatched = p.Unmatched;
            
            p = inputParser;
            addParameter(p, 'Covariance', speye(obj.NumObs));
            
            parse(p, unMatched);
            s = p.Results;
            obj.Covariance = s.Covariance;
        end

        function y = observeWithError(obj, t, x)
            y = obj.ErrorModel.addError(t, obj.observeWithoutError(t, x));
        end

        function y = observeWithoutError(~, ~, x)
            y = x;
        end

        function H = linearization(obj, ~)
            H = speye(obj.NumVars);
        end

    end

end
