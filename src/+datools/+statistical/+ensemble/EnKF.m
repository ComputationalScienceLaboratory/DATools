classdef EnKF < handle
    
    properties
        Model
        ModelError
        Observation
        Ensemble
        Inflation
        Localization
        Parallel
        RIPIterations
    end
    
    properties (Dependent)
        BestEstimate
        NumEnsemble
    end
    
    properties (Hidden = true)
        DidRIP
    end
    
    methods
        function obj = EnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            addParameter(p, 'ModelError', datools.error.Error);
            addParameter(p, 'NumEnsemble', 1);
            addParameter(p, 'Inflation', 1);
            addParameter(p, 'Localization', @(~, ~, H) datools.tapering.trivial(H));
            addParameter(p, 'Parallel', false);
            addParameter(p, 'RIPIterations', 0);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.Model        = s.Model;
            obj.ModelError   = s.ModelError;
            obj.Inflation    = s.Inflation;
            obj.Localization = s.Localization;
            obj.Parallel     = s.Parallel;
            obj.RIPIterations = s.RIPIterations;
            obj.DidRIP        = false;
            ensN = s.NumEnsemble;
            
            kept = p.Unmatched;
            
            p = inputParser;
            addParameter(p, 'Observation', datools.observation.Observation(s.Model.NumVars));
            addParameter(p, 'EnsembleGenerator', @(x) randn(s.Model.NumVars, x));
            parse(p, kept);
            
            s = p.Results;
            
            obj.Ensemble = s.EnsembleGenerator(ensN);
            obj.Observation = s.Observation;
            
        end
        
        function forecast(obj)
            
            times = zeros(obj.NumEnsemble, 1);
            
            if obj.Parallel
                rhs = obj.Model.ODEModel.Rhs.F;
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
            
            if obj.DidRIP
                ripits = 0;
            else
                ripits = obj.RIPIterations;
            end
            
            for ripit = 1:(ripits + 1)
                
                inflation = obj.Inflation;
                
                tc = obj.Model.TimeSpan(1);
                
                xf = obj.Ensemble;
                ensN = obj.NumEnsemble;
                
                
                xfm = mean(xf, 2);
                
                Af = xf - repmat(xfm, 1, ensN);
                Af = inflation*Af;
                
                xf = repmat(xfm, 1, ensN) + Af;
                
                Hxf = obj.Observation.observeWithoutError(tc, xf);
                Hxfm = mean(Hxf, 2);
                

                
                HAf = Hxf - repmat(Hxfm, 1, ensN);
                
                % tapering
                
               
                
                if isempty(obj.Localization)
                    rhoHt  = ones(size(Af, 1), size(HAf, 1));
                    HrhoHt = ones(size(HAf, 1), size(HAf, 1));
                else
                    H = obj.Observation.linearization(tc, xfm);
                    rhoHt = obj.Localization(tc, xfm, H);
                    HrhoHt = H*rhoHt;
                end
                
                
                PfHt = rhoHt.*((1/(ensN - 1))*(Af*(HAf')));
                HPfHt = HrhoHt.*((1/(ensN - 1))*(HAf*(HAf')));
                
                S = HPfHt + R;
                
                dS = decomposition(S, 'lu');
                
                d = y - Hxfm;
                
                xam = xfm + PfHt*(dS\d);
                
                Aa = Af + PfHt*(dS\(sqrtm(R)*randn(size(HAf)) - HAf));
                
                %Aa = inflation*Aa;
                
                obj.Ensemble = repmat(xam, 1, ensN) + Aa;
                
                obj.Model.update(0, obj.BestEstimate);
                
            end
            
            obj.DidRIP = true;
            
        end
        
        function ensN = get.NumEnsemble(obj)
            ensN = size(obj.Ensemble, 2);
        end
        
        
        function ensN = get.BestEstimate(obj)
            ensN = mean(obj.Ensemble, 2);
        end
        
        function setMean(obj, xam)
            
            Xf  = obj.Ensemble;
            
            ensN = size(Xf, 2);
            
            xfm  = mean(Xf,  2);
            
            Af  = (Xf  - repmat(xfm,  1, ensN));
            
            Xf  = repmat(xam, 1, ensN) + Af;
            
            obj.Ensemble= Xf;
            
            
        end
        
        function scaleAnomalies(obj, scale)
            
            Xf  = obj.Ensemble;
            
            ensN = size(Xf, 2);
            
            xfm  = mean(Xf,  2);
            
            Af  = scale*(Xf  - repmat(xfm,  1, ensN));
            
            Xf  = repmat(xfm, 1, ensN) + Af;
            
            obj.Ensemble= Xf;
            
            
        end
        
    end
    
end