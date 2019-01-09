classdef DEnKF < handle
    
    properties
        Model
        ModelError
        Observation
        Ensemble
        Inflation
        Localization
        Parallel
    end
    
    properties (Dependent)
        BestEstimate
        NumEnsemble
    end
    
    methods
        function obj = DEnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'csl.datools.Model'));
            addParameter(p, 'ModelError', csl.datools.error.Error);
            addParameter(p, 'NumEnsemble', 1);
            addParameter(p, 'Inflation', 1);
            addParameter(p, 'Localization', @(~, ~, H) csl.datools.tapering.trivial(H));
            addParameter(p, 'Parallel', false);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.Model        = s.Model;
            obj.ModelError   = s.ModelError;
            obj.Inflation    = s.Inflation;
            obj.Localization = s.Localization;
            obj.Parallel     = s.Parallel;
            ensN = s.NumEnsemble;
            
            kept = p.Unmatched;
            
            p = inputParser;
            addParameter(p, 'Observation', csl.datools.observation.Observation(s.Model.NumVars));
            addParameter(p, 'EnsembleGenerator', @(x) randn(s.Model.NumVars, x));
            parse(p, kept);
            
            s = p.Results;
            
            obj.Ensemble = s.EnsembleGenerator(ensN);
            obj.Observation = s.Observation;
            
        end
        
        function forecast(obj)
            
            times = zeros(obj.NumEnsemble, 1);
            
            if obj.Parallel
                rhs = obj.Model.ODEModel.F;
                tspan = obj.Model.TimeSpan;
                solver = obj.Model.Solver;
                ens    = obj.Ensemble;
                ensN   = obj.NumEnsemble;
                
                parfor ensi = 1:ensN
                    
                    [t, y] = solver(rhs, tspan, ens(:, ensi));
                    
                    time = t(end) - t(1);
                    yend = y(end, :).';
                    
                    ens(:, ensi) = yend;
                    times(ensi) = time;
                end
                
                for ensi = 1:ensN
                    obj.Ensemble(:, ensi) = obj.ModelError.adderr(obj.Model.TimeSpan(end), ens(:, ensi));
                end
                
            else
                for ensi = 1:obj.NumEnsemble
                    [time, yend] = obj.Model.solve([], obj.Ensemble(:, ensi));
                    
                    obj.Ensemble(:, ensi) = obj.ModelError.adderr(obj.Model.TimeSpan(end), yend);
                    times(ensi) = time;
                end
            end
            
            obj.Model.update(mean(times), obj.BestEstimate);
            
        end
        
        function analysis(obj, R, y)
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            
            xfm = mean(xf, 2);
            
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            
            Af = xf - repmat(xfm, 1, ensN);
            
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            
            % tapering
            inflation = obj.Inflation;
            H = obj.Observation.linearization(tc, xf);
            rhoHt = obj.Localization(tc, xfm, H);
            HrhoHt = H*rhoHt;
            
            
            PfHt = rhoHt.*((1/(ensN - 1))*(Af*(HAf')));
            HPfHt = HrhoHt.*((1/(ensN - 1))*(HAf*(HAf')));
            
            S = HPfHt + R;
            
            d = y - Hxfm;
            
            xam = xfm + PfHt*(S\d);
            
            Aa = Af - 0.5*PfHt*(S\HAf);
            
            Aa = inflation*Aa;
            
            obj.Ensemble = repmat(xam, 1, ensN) + Aa;
            
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
        function ensN = get.NumEnsemble(obj)
            ensN = size(obj.Ensemble, 2);
        end
        
        
        function ensN = get.BestEstimate(obj)
            ensN = mean(obj.Ensemble, 2);
        end
        
    end
    
end