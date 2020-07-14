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

            % abuse 
            inflation = obj.Inflation;
            M = obj.SurrogateEnsN;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            AeqT = [kron(speye(ensN), ones(1, ensN + M)); kron(ones(ensN, 1), speye(ensN + M)).'];
            lbT = zeros((ensN + M)*ensN, 1);
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
            
            if M == 0
                gamma = 0;
            end
            mu = trC/n;
            
            %mu =  sum(s1.^2)/sum(s2.^2);
            Asynth = sqrt(mu)*Bsqrtf*randn(n, M);
            
            laplace = obj.Laplace;
            if laplace
                Z = exprnd(1, 1, M);
                Asynth = sqrt(Z).*Asynth;
            end
            
            Asynth = Asynth - mean(Asynth, 2);
            
            Afrak = [Af, inflation*Asynth];
            chiF = xfm + Afrak;
            
            Hchif = obj.Observation.observeWithoutError(tc, chiF);
            
            xdist = zeros(size(chiF));
            w = zeros(ensN + M, 1);
            
            for i = 1:(ensN + M)
                X_temp = chiF - chiF(:, i);
                xdist(i, :) = vecnorm(X_temp, 2).^2;
                inn = y - Hchif(:, i);
                w(i) = exp(-0.5*inn'*(R\inn));
                
            end
            
            xdist = xdist(:, 1:ensN);
            
            w(1:ensN) = w(1:ensN)*(1 - gamma)/(ensN);
            w((ensN + 1):end) = w((ensN + 1):end)*gamma/(M);
            
            w = w/sum(w);
            
            beqT = [ones(ensN, 1)/ensN; w];
            f = xdist(:);
            Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
            Tx = ensN*reshape(Tx, ensN + M, ensN);
            
            xa = chiF*(Tx);
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end
