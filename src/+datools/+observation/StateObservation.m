classdef StateObservation < handle
    %STATEOBSERVATION this is the class that handles the observed state
    %   Observed state are noisy snapshots of reality (truth)

    properties
        Y         % Current Observation of States (maybe sparse)
        R         % Observation Error Covariance
        Noise     % Noise
        NumObs    % Number of Observations(Make sure > 0)
    end
    
    methods
        function obj = StateObservation(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Y', []);
            addParameter(p, 'Noise', datools.error.Error);
            validationFunction = @(x) (x>0);
            addParameter(p, 'NumObs', 1, validationFunction);
            
            parse(p, varargin{:});

            s = p.Results;
            
            
            obj.Y = s.Y;
            obj.Noise = s.Noise;
            obj.NumObs = s.NumObs;
            
            unMatched = p.Unmatched;
            
            p = inputParser;
            addParameter(p, 'R', speye(obj.NumObs));
            
            parse(p, unMatched);
            s = p.Results;
            obj.R = s.R;
        end
        
        function updateY(obj, y)
            obj.Y = y;
        end
        
        function addNoise(obj, t, ~)
            obj.Y = obj.Noise.adderr(t, obj.Y);
        end
    end

end