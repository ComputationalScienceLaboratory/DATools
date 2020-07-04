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
            obj.Bsqrt = sqrtm(s.B);
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
            
            Covsqrt = (1/sqrt(ensN - 1)*Af);           
            Bsqrtf = obj.Bsqrt;
     
            s = svd(Bsqrtf\Covsqrt);
            trC = sum(s.^2);
            tr2C = trC*trC;
            trC2 = sum(s.^4);
            
            n = size(xf, 1);
            
            
            NN = ensN - 1;
            gamma = min(1, ((NN - 2)/NN * trC2 + tr2C)/((NN + 2)*(trC2 - tr2C/n)));
            mu = trC/n;
            
            %mu =  sum(s1.^2)/sum(s2.^2);
            Asynth = sqrt(mu)*Bsqrtf*randn(n, M);
            
            laplace = obj.Laplace;
            if laplace
                Z = exprnd(1, 1, M);
                Asynth = Z.*Asynth;
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
