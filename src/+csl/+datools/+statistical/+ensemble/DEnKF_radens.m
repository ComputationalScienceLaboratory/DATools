classdef DEnKF_radens < handle
    
    properties
        Model
        ModelError
        Observation
        Ensemble
        Inflation
        Localization
        Parallel
        RadEns
        Step
    end
    
    properties (Dependent)
        BestEstimate
        NumEnsemble
    end
    
    methods
        function obj = DEnKF_radens(varargin)
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
            
            radsamples = 5;
            radens = 10;
            %obj.RadEns = sqrt(repmat(linspace(0.15, 30, radens), radsamples, 1)).*randn(radsamples, radens);
            
            %obj.RadEns = sqrt(repmat(linspace(0.15, 200, radens), radsamples, 1));
            
            obj.RadEns = 4*randn(radsamples, radens);
            
            obj.Step = 0;
            
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
            
            d = y - Hxfm;
            
            Af = xf - repmat(xfm, 1, ensN);
            
            HAf = Hxf - repmat(Hxfm, 1, ensN);
            
            inflation = obj.Inflation;
            H = obj.Observation.linearization(tc, xf);
            
            rf = obj.RadEns;
            
            sfm = mean(rf, 2);
            
            rin = 1.0;
            
            Sf = rf - repmat(sfm, 1, size(rf, 2));
            
            Sf = Sf;
            
            %SPf = (Sf*Sf')/(size(rf, 2) - 1);
                        
            xams = [];
            
            rfs = mean(rf.^2, 1);
            
            if obj.Parallel
                nworkers = 0;
            else
                nworkers = 0;
            end
            
            parfor (i = 1:numel(rfs), nworkers)
                
                % tapering
                rhoHt = obj.Localization(tc, xfm, H, rfs(i));
                HrhoHt = H*rhoHt;
                
                
                PfHt = rhoHt.*((1/(ensN - 1))*(Af*(HAf')));
                HPfHt = HrhoHt.*((1/(ensN - 1))*(HAf*(HAf')));
                
                S = HPfHt + R;
                
                xam = xfm + PfHt*(S\d);
                
                xams = [xams, xam];

            end
            
            xamm = mean(xams, 2);
            
            Aam = xams - repmat(xamm, 1, size(rf, 2));
            
            
            
            SPfHt = rin*(1/(ensN - 1))*(Sf*((H*Aam)'));
            HSPfHt = (1/(ensN - 1))*((H*Aam)*((H*Aam)'));
            
            dS = y - H*xamm;
            
            SS = HSPfHt + R;
            
            sam = sfm + SPfHt*(SS\dS);
            
            Sa = Sf - 0.5*SPfHt*(SS\(H*Aam));
            
            obj.RadEns = repmat(sam, 1, size(rf, 2)) + Sa;
            
            r = mean(mean(obj.RadEns.^2, 1));
            
            rvar = var(mean(obj.RadEns.^2, 1));
            
            fprintf('r: %.5f,       var: %.5f     ', r, rvar);
            
            %pause(1)
            
            rhoHt = obj.Localization(tc, xfm, H, r);
            HrhoHt = H*rhoHt;
            PfHt = rhoHt.*((1/(ensN - 1))*(Af*(HAf')));
            HPfHt = HrhoHt.*((1/(ensN - 1))*(HAf*(HAf')));  
            S = HPfHt + R;
            xam = xfm + PfHt*(S\d);
            
            Aa = Af - 0.5*PfHt*(S\HAf);
            
            Aa = inflation*Aa;
            
            obj.Step = obj.Step + 1;
            
            % add process noise to radius
            
            if obj.Step < 500
                mu = 0.25;
            else
                mu = 0.01;
            end
            %Qsqrt = sqrt(mu)*eye(size(obj.RadEns, 1));
            
            if mu > 0
                Q = mu*spdiags(repmat([-1/2, 1, -1/2], size(obj.RadEns, 1), 1), [-1, 0, 1], size(obj.RadEns, 1), size(obj.RadEns, 1));
                Qsqrt = sqrtm(full(Q));
                
                
                obj.RadEns = obj.RadEns + mu*Qsqrt*randn(size(obj.RadEns));
                
            end
            
            
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