classdef MLEnKF < handle
    
    properties
        Models
        ModelError
        Observation
        Ensembles
        Inflation
        Localization
        Parallel
    end
    
    properties (Dependent)
        BestEstimate
        NumEnsemble
    end
    
    methods
        function obj = MLEnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Models', @(x) isa(x, 'cell'));
            addParameter(p, 'ModelError', csl.datools.error.Error);
            addParameter(p, 'NumEnsemble', 1);
            addParameter(p, 'Inflation', 1);
            addParameter(p, 'Localization', @(~, ~, H) csl.datools.tapering.trivial(H));
            addParameter(p, 'Parallel', false);
            parse(p, varargin{:});
            
            s = p.Results;
            
            obj.Models       = s.Models;
            obj.ModelError   = s.ModelError;
            obj.Inflation    = s.Inflation;
            obj.Localization = s.Localization;
            obj.Parallel     = s.Parallel;
            ensN = s.NumEnsemble;
            
            kept = p.Unmatched;
            
            p = inputParser;
            addParameter(p, 'Observation', csl.datools.observation.Observation(s.Models{1}.NumVars));
            addParameter(p, 'EnsembleGenerator', @(x) randn(s.Models{1}.NumVars, x));
            parse(p, kept);
            
            s = p.Results;
                        
            for ei = 1:(numel(obj.Models)*2 - 1)
                obj.Ensembles{ei} = s.EnsembleGenerator(ensN);
            end
            
            obj.Observation = s.Observation;
            
        end
        
        function forecast(obj)
            
            
            % We will start from the lowest level
            for ei = 1:numel(obj.Ensembles)
                mi = floor((ei + 1)/2);
                
                times = zeros(obj.NumEnsemble, 1);
                
                if obj.Parallel
                    rhs = obj.Models{mi}.ODEModel.F;
                    tspan = obj.Models{mi}.TimeSpan;
                    solver = obj.Models{mi}.Solver;
                    %ens    = obj.Ensembles{ei};
                    % test
                    if ei == 2
                        ens    = obj.Ensembles{3};
                    else
                        ens    = obj.Ensembles{ei};
                    end
                    ensN   = obj.NumEnsemble;
                    
                    parfor ensi = 1:ensN
                        
                        [t, y] = solver(rhs, tspan, ens(:, ensi));
                        
                        time = t(end) - t(1);
                        yend = y(end, :).';
                        
                        ens(:, ensi) = yend;
                        times(ensi) = time;
                    end
                    
                    for ensi = 1:ensN
                        obj.Ensembles{ei}(:, ensi) = obj.ModelError.adderr(obj.Models{mi}.TimeSpan(end), ens(:, ensi));
                    end
                    
                else
                    for ensi = 1:obj.NumEnsemble
                        [time, yend] = obj.Models{mi}.solve([], obj.Ensembles{ei}(:, ensi));
                        
                        obj.Ensembles{ei}(:, ensi) = obj.ModelError.adderr(obj.Models{mi}.TimeSpan(end), yend);
                        times(ensi) = time;
                    end
                end
                
                if mod(ei, 2) == 1
                    obj.Models{mi}.update(mean(times), obj.BestEstimate);
                end
                
            end
            
        end
        
        function analysis(obj, R, y)
            
            tc = obj.Models{1}.TimeSpan(1);
            
            H = obj.Observation.linearization(tc, []);
            
            
            xfm  = cell(numel(obj.Ensembles), 1);
            Hxfm = cell(numel(obj.Ensembles), 1);
            Af   = cell(numel(obj.Ensembles), 1);
            HAf  = cell(numel(obj.Ensembles), 1);
            
            PfHt = 0;
            HPfHt = 0;
            
            if tc < 100
                ssmall = 1/2;
            else
                ssmall = 1/2;
            end
            
            for ei = 1:numel(obj.Ensembles)
                
                mi = floor((ei + 1)/2);
                
                xf = obj.Ensembles{ei};
                ensN = obj.NumEnsemble;
                
                xfm{ei} = mean(xf, 2);
                
                Hxf = obj.Observation.observeWithoutError(tc, xf);
                Hxfm{ei} = mean(Hxf, 2);
                
                Af{ei} = xf - repmat(xfm{ei}, 1, ensN);
                
                HAf{ei} = Hxf - repmat(Hxfm{ei}, 1, ensN);
                
                s = mod(ei, 2) * 2 - 1;
                
                if ei ~= numel(obj.Ensembles)
                    s = s*ssmall;
                end
                

                
                rhoHt = obj.Localization{mi}(tc, xfm{ei}, H);
                HrhoHt = H*rhoHt;
                
                PfHt = PfHt + s*rhoHt.*((1/(ensN - 1))*(Af{ei}*(HAf{ei}.')));
                tmp = ((1/(ensN - 1))*(HAf{ei}*(HAf{ei}.')));
                %if ei ~= numel(obj.Ensembles)
                %    tmp = tmp - diag(diag(tmp));
                %end
                HPfHt = HPfHt + s*HrhoHt.*tmp;
                
            end
            
            inflation = obj.Inflation;
            
            % tapering
            
            %rhoHt = obj.Localization(tc, xfm{ei}, H);
            %HrhoHt = H*rhoHt;
            
            % do Kody Law's thing
            %[Ud, Sd, Vd] = svd(HPfHt);
            %k = sum(diag(Sd) > 0);
            %HPfHt = Ud(:, 1:k)*Sd(1:k, 1:k)*(Vd(:, 1:k).');
            
            %PfHt  = rhoHt.*PfHt;
            %HPfHt = HrhoHt.*HPfHt;
            
            dHPfHt = diag(HPfHt);
            HPfHt = HPfHt - diag(dHPfHt);
            dHPfHt(dHPfHt < 0) = 0;
            HPfHt = HPfHt + diag(dHPfHt);
            
            %din = min(min(diag(HPfHt)), 0);
            %HPfHt = HPfHt - din*eye(size(HPfHt));
            
            
            S = HPfHt + R;
            
            %lambda = 1;
            %
            %             [~, flag] = chol(HPfHt + R);
            %
            %             while flag
            %                 [~, flag] = chol(HPfHt + lambda*R);
            %                 lambda = lambda*1.25;
            %             end
            
            %S = HPfHt + lambda*R;
            
            yp  = cell(numel(obj.Models), 1);
            
            for mi = 1:numel(obj.Models)
                yp{mi} = sqrtm(R)*randn(size(HAf{1}));
            end
            
            for ei = 1:numel(obj.Ensembles)                
                d = y - Hxfm{ei};
                
                mi = floor((ei + 1)/2);
                
                xam = xfm{ei} + PfHt*(S\d);
                
                Aa = Af{ei} + PfHt*(S\(yp{mi} - HAf{ei}));
                
                Aa = inflation*Aa;
                
                obj.Ensembles{ei} = repmat(xam, 1, ensN) + Aa;
                
                obj.Models{mi}.update(0, obj.BestEstimate);
                
            end
            
            
        end
        
        function ensN = get.NumEnsemble(obj)
            ensN = size(obj.Ensembles{1}, 2);
        end
        
        
        function xa = get.BestEstimate(obj)
            
            xa = 0;
            
            for ei = 1:numel(obj.Ensembles)
                
                s = mod(ei, 2) * 2 - 1;
                
                xa = xa + s*mean(obj.Ensembles{ei}, 2);
                
            end
            
            
        end
        
    end
    
end