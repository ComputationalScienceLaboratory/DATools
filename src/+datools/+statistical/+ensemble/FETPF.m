classdef FETPF < datools.statistical.ensemble.EnF
    
    properties
        B
        Bsqrt
        SurrogateEnsN
        Laplace
    end
    
    methods
        function obj = FETPF(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'B', 1);
            addParameter(p, 'SurrogateEnsembleSize', 2);
            addParameter(p, 'Laplace', false);
            
            parse(p, varargin{2:end});
            
            s = p.Results;
            
            kept = p.Unmatched;
            
            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);
            
            obj.B = s.B;
            if iscell(obj.B)
                Bs = cell(numel(obj.B), 1);
                for i = 1:numel(obj.B)
                    Bs{i} = sqrtm(obj.B{i});
                end
                obj.Bsqrt = Bs;
            else
                obj.Bsqrt = sqrtm(s.B);
            end
            obj.SurrogateEnsN = s.SurrogateEnsembleSize;
            obj.Laplace = s.Laplace;
            
        end
    end
    
    methods
        
        function analysis(obj, R, y)
            
            inflation = obj.Inflation;
            tau = obj.Rejuvenation;
            
            ensM = obj.SurrogateEnsN;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            dR = decomposition(R, 'chol');
            
            AeqT = [kron(speye(ensN), ones(1, ensN + ensM)); kron(ones(ensN, 1), speye(ensN + ensM)).'];
            lbT = zeros((ensN + ensM)*ensN, 1);
            optsT = optimoptions('linprog', 'Display', 'off');
            
            
            xfm = mean(xf, 2);
            Af = xf - xfm;
            
            n = size(xf, 1);
            
            Covsqrt = (1/sqrt(ensN - 1)*Af);
            Bsqrtf_all = obj.Bsqrt;
            
            NN = ensN - 1;
            
            if iscell(Bsqrtf_all)
                
                Bi = 0;
                Uhopt = inf;
                
                for i = 1:numel(Bsqrtf_all)
                    
                    Bsqrtfi = Bsqrtf_all{i};
                    s = svd(Bsqrtfi\Covsqrt);
                    trC = sum(s.^2);
                    tr2C = trC*trC;
                    trC2 = sum(s.^4);
                    
                    % calculate sphericity
                    Uh = (n*trC2/tr2C - 1)/(n - 1);
                    
                    if Uh < Uhopt
                        Uhopt = Uh;
                        Bi = i;
                    end
                    
                end
                Bsqrtf = Bsqrtf_all{Bi};
                
            else
                
                Bsqrtf = Bsqrtf_all;
                s = svd(Bsqrtf\Covsqrt);
                trC = sum(s.^2);
                tr2C = trC*trC;
                trC2 = sum(s.^4);
                
                % calculate sphericity
                Uhopt = (n*trC2/tr2C - 1)/(n - 1);
                
            end
            
            gamma = min((NN - 2)/(NN*(NN + 2)) + ((n + 1)*NN - 2)/(Uhopt*NN*(NN + 2)*(n - 1)), 1);
            
            if ensM == 0
                gamma = 0;
            end
            mu = trC/n;
            
            %mu =  sum(s1.^2)/sum(s2.^2);
            Asynth = sqrt(mu)*Bsqrtf*randn(n, ensM);
            
            laplace = obj.Laplace;
            if laplace
                Z = exprnd(1, 1, ensM);
                Asynth = sqrt(Z).*Asynth;
            end
            
            Asynth = Asynth - mean(Asynth, 2);
            
            Afrak = [Af, inflation*Asynth];
            chiF = xfm + Afrak;
            
            Hchif = obj.Observation.observeWithoutError(tc, chiF);
            
            xdist = zeros(size(chiF));
            
            for i = 1:(ensN + ensM)
                X_temp = chiF - chiF(:, i);
                xdist(i, :) = vecnorm(X_temp, 2).^2;
            end
            
            t0 = Hchif - y;
            
            % more efficient way of calculating weights
            as = (-0.5*sum(t0.*(dR\t0), 1)).';
            m = max(as);
            w = exp(as - (m + log(sum(exp(as - m)))));
            
            xdist = xdist(:, 1:ensN);
            
            w(1:ensN) = w(1:ensN)*(1 - gamma)/(ensN);
            w((ensN + 1):end) = w((ensN + 1):end)*gamma/(ensM);
            
            w = w/sum(w);
            
            beqT = [ones(ensN, 1)/ensN; w];
            f = xdist(:);
            Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
            Tx = ensN*reshape(Tx, ensN + ensM, ensN);
            
            xa = chiF*(Tx);
            
            obj.Ensemble = xa;
            obj.Weights = ones(ensN, 1)/ensN;
            
            % rejuvenation of this sort is not performed by default,
            % however it is included in this algorithm purely for
            % completeness
            obj.rejuvenate(tau);
            
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
