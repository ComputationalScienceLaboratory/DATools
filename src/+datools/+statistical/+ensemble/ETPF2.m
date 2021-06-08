classdef ETPF2 < datools.statistical.ensemble.EnF
    
    methods
        
        function analysis(obj, R, y)

            tau = obj.Rejuvenation;
            
            tc = obj.Model.TimeSpan(1);
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            AeqT = [kron(speye(ensN), ones(1, ensN)); kron(ones(ensN, 1), speye(ensN)).'];
            lbT = zeros((ensN)*ensN, 1);
            optsT = optimoptions('linprog', 'Display', 'off');
           
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            
            xdist = zeros(ensN, ensN);
            
            for i = 1:ensN
                xtemp = xf - repmat(xf(:, i), 1, ensN);
                xdist(i, :) = vecnorm(xtemp).^2;
            end
            
            t0 = Hxf - y;
            w = exp(-0.5*sum(t0.*(R\t0), 1));
            w = w/sum(w);
            
            beqT = [ones(ensN, 1)/ensN; w];
            f = xdist(:);
            Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
            Tx = ensN*reshape(Tx, ensN, ensN);
            
            
            % Delta
            W = diag(w);
            Bric = Tx - w*ones(1, ensN);
            Aric = ensN*(W - w*w.') - Bric*Bric.';
            
            rs = @(D) reshape(D, ensN, ensN);
            urs = @(D) reshape(D, ensN*ensN, 1);
            f = @(t, D) urs(-Bric * rs(D) - rs(D)*(Bric.') + Aric - rs(D)*rs(D));
            D = zeros(ensN, ensN);
            D = rs(dp(f, D(:)));
            
            
            P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
            xa = xf*(Tx + D + P);
            
            obj.Ensemble = xa;
            obj.Model.update(0, obj.BestEstimate);
            
        end
        
    end
    
end


function y = dp(f, y0)

n = sqrt(numel(y0));

h = 0.1;
y = y0;

abstol = 1e-6;
reltol = 1e-6;

ks = zeros([numel(y0), 7]);

A = [0, 0, 0, 0, 0, 0, 0; ...
    1/5, 0, 0, 0, 0, 0, 0; ...
    3/40, 9/40, 0, 0, 0, 0, 0; ...
    44/45, -56/15, 32/9, 0, 0, 0, 0; ...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0; ...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0; ...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];

bs = A(end, :);
bhs = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];

t = 0;

while true    
    for s = 1:7
        ks(:, s) = h*f(0, y + sum(ks(:, 1:(s - 1)).*A(s, 1:(s - 1)), 2));
    end
    
    
    
    ynew = y + sum(ks.*bs, 2);
    yhnew = y + sum(ks.*bhs, 2);
    
    sc = abstol + max(abs(ynew), abs(yhnew))*reltol;
    err = rms((ynew - yhnew)./sc);
    
    orderE = 4;
    fac = 0.38^(1/(orderE + 1));
    
    facmin = 0.2;
    
    if err > 1 || isnan(err)
        facmax = 1;
        
        hnew = h*min(facmax, max(facmin, fac*(1/err)^(1/(orderE + 1))));
    else
               
        % condition used by Acevedo, de Wiljes and Reich.
        if norm(reshape(ynew, n, n) - reshape(y, n, n), 'inf') < 1e-3 || t > 1000
            break;
        end
        
        y = ynew;
        t = t + h;
        
        
        facmax = 3;
        
        hnew = h*min(facmax, max(facmin, fac*(1/err)^(1/(orderE + 1))));

    end

        
    h = hnew;
    
    
end

end

