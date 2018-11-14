function pdawnet
clc;
% for L = [1 4 10 30]
%     for N = [1 2 5 10]
%         for Nchoose = [1 2 5 10 50]

zzz = load('lorenz96bigmodelobsforenet2.mat');
net = zzz.net;
for L = [8]
    for N = [10]
        for Nchoose = [10]
            
            fprintf('L: %d, N: %d, Nc: %d\n', L, N, Nchoose);
            
            rng(6);
            
            nvar = 40;
            nvarnature = 2*nvar;
            
            %nvarnature = nvar;
            
            if nvarnature == nvar
                naturetotruth = speye(nvar);
            else
                naturetotruth = spdiags(repmat([1/4, 1/4, 1/2, 1/4, 1/4], nvarnature, 1), [-nvarnature+1, -1, 0, 1, nvarnature-1], nvarnature, nvarnature);
                naturetotruth = naturetotruth(1:2:nvarnature, :);
            end

            %spy(naturetotruth)
            
            %modeltruth = csl.odetestproblems.lorenz96.presets.Canonical;
            modeltruth = csl.odetestproblems.lorenz96.presets.RandomIC(nvarnature, 6);
            modelmodel = csl.odetestproblems.lorenz96.presets.RandomIC(nvar, 6);
            
            modelmodel.Y0 = 4*randn(modelmodel.NumVars, 1);
            
            Esqrt = sqrtm(speye(modelmodel.NumVars));
            
            dt = 0.05;
            
            %ti = @(f, t, y) csl.odeutils.rk4(f, t, y, dt);
            
            ti = @(f, t, y) rk4(f, t, y, dt);
            
            ntol = 1e-8;
            nJvp = modeltruth.JacobianVectorProduct;
            solvernature = @(f, t, y) csl.utils.RODAS4(f, t, y, ntol, nJvp);

            
            %dt = 0.01;
            %ti = @(f, t, y) fe(f, t, y, 1);
            
            % this is an estimate
            %Q = 2*speye(nvar);
            Q = 0*speye(nvar);
            
            Qsqrt = sqrtm(Q);
            
            sdallow = 4;
            
            maxits = 1024;
            
            % at which point to calculate the analysis error
            cea = floor((L + 1)/2);
            
            
            xf = repmat(modelmodel.Y0, 1, N) + Esqrt*randn(nvar, N);
            
            % propogate the trajectories
            
            for j = 1:20
                [~, Y] = ti(modeltruth.F, [0 dt], modeltruth.Y0);
                modeltruth.Y0 = Y(end, :).';
                
                for ensi = 1:N
                    [~, Y] = ti(modelmodel.F, [0 dt], xf(:, ensi));
                    xf(:, ensi) = Y(end, :).';
                end
                modelmodel.Y0 = mean(xf, 2);
            end
            
            obsvars = [2:2:(nvar/2 - 1), nvar/2:nvar];
            %obsvars = 1:nvar;
            I = speye(modelmodel.NumVars);
            H = I(obsvars, :);
            
            
            spinup = 10;
            time = spinup*11;
            
            err = zeros(time, 1);
            
            y = [];
            xt = [];
            
            R = speye(numel(obsvars));
            Rsqrt = sqrtm(R);
            
            dispstr = '';
            
            xfprev = [];
            
            for i = 1:time
                
                
                yprev = y;
                
                % if we don't have enough points in time, we do not discard, else we do
                % The size of y is n by L+1, as L is the over sampling factor
                if size(y, 2) > L
                    y = y(:, 2:end);
                    xt = xt(:, 2:end);
                end
                
                %errf = rms(modelmodel.Y0 - modeltruth.Y0)
                
                xt = [xt, naturetotruth*modeltruth.Y0];
                % observe the full trajectory with error
                y = [y, (H*(xt(:, end)) + Rsqrt*randn(numel(obsvars), 1))];
                
                % perform PDA, hugs and kisses.
                
                if size(y, 2) > L && i > L + 1
                    
                    
                    %opts = optimoptions('fminunc', 'Display', 'off', 'MaxIterations', inf);
                    
                    opts = optimoptions('fminunc', 'Display', 'off', 'algorithm', 'quasi-newton', 'MaxIterations', 10);
                    
                    
                    %opts = optimoptions('fminunc', 'Display', 'off', ...
                    %    'algorithm', 'quasi-newton', 'MaxIterations', 1024, 'HessUpdate', 'bfgs', 'SpecifyObjectiveGradient', true);
                    %opts = optimoptions('fminunc', 'Display', 'off', 'algorithm', 'quasi-newton', ...
                    %    'MaxIterations', inf, 'HessUpdate', 'bfgs', 'SpecifyObjectiveGradient', true);
                    
                    
                    %opts = optimoptions('fmincon', 'Display', 'off', ...
                    %    'MaxIterations', 256, 'HessianApproximation', {'lbfgs', max([L, 20])}, 'SpecifyObjectiveGradient', true);
                    %opts = optimoptions('fmincon', 'Display', 'off', ...
                    %    'MaxIterations', inf, 'HessianApproximation', 'lbfgs', 'SpecifyObjectiveGradient', false);
                    
                    
                    %opts = optimoptions('fminunc', 'Display', 'off', 'algorithm', 'quasi-newton', 'MaxIterations', inf, 'HessUpdate', 'bfgs');
                    
                    %opts = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', inf);
                    
                    %opts = optimoptions('patternsearch', 'Display', 'off', 'MaxIterations', 100);
                    
                    %xens = zeros(modelmodel.NumVars, N);
                    
                    nv = modelmodel.NumVars;
                    
                    likelihoods = zeros(N, 1);
                    
                    
                    F = modelmodel.F;
                    Javp = modelmodel.JacobianAdjointVectorProduct;
                    
                    
                    parfor ensi = 1:N
                        xens = reshape(xf(:, ensi, :), nv, L + 1);
                        %yens = (y + Rsqrt*randn(numel(obsvars), L + 1));
                        yens = y;
                        
                        qens = Qsqrt*randn(nv, L + 1);
                        
                        cf = @(yy) cost(yy, qens, L, ti, dt, F, Javp);
                        
                        sc = @(yy) stop(yy, y, H, R, L, sdallow);
                        
                        
                        % We predict the current best trajectory using a neural
                        % net
                        
                        yprevens = yprev;
                        
                        xfprevens = reshape(xfprev(:, ensi, :), nv, L + 1);
                        
                        %netin = [yprevens - H*xfprevens; xfprevens];
                        netin = [xfprevens; yprevens - H*xfprevens];
                        
                        d = net(netin);
                        
                        % weight the observations by likelyhood;
%                         for l = 1:(L + 1)
%                             lhe = lhe + 0.5 * (H*xens(:, l) - yensc(:, l))'*(R\(H*xens(:, l) - yensc(:, l)));
%                             lhn = lhn + 0.5 * (H*xens(:, l) - yensnet(:, l))'*(R\(H*xens(:, l) - yensnet(:, l)));
%                         end
                        
                        lhe = 0;
                        lhn = 1;
                        
                        %yens = (lhe*yensc + lhn*yensnet)/(lhe + lhn);
                        
                        
                        
                        %d = yens - H*xens;
                        
                        xenswo = xens + H\(d);
 
                        
                        
%                         [~, grad] = cf(xenswo(:));
%                         fdgrad = zeros(size(grad));
%                         for jj = 1:numel(fdgrad)
%                             e = zeros(size(grad));
%                             e(jj) = 1;
%                             h = sqrt(eps);
%                             fdgrad(jj) = (cf(xenswo(:) + h*e) - cf(xenswo(:) - h*e))/(2*h);
%                         end
%                         norm(grad - fdgrad)/norm(fdgrad)
%                         return;
                        %xa = fminunc(cf, xenswo(:), opts);
                        
                        %xa = fmincon(cf, xenswo(:),[],[],[],[],[],[],[], opts);
                        
                        xa = csl.utils.gd(cf, xenswo(:), maxits, eps, eps, sc);
                        
                        %xa = patternsearch(cf, xenswo(:),[],[],[],[],[],[],[], opts);
                        
                        xa = reshape(xa, nv, 1, L + 1);
                        %xa = xa(:, 1, :);
                        
                        xf(:, ensi, :) = xa;
                        %xf(:, ensi, L + 1) = xa(:, :, L + 1);
                        
                        
                        xfc = reshape(xf(:, ensi, :), nv, L + 1);
                        
                        for l = 1:(L + 1)
                            likelihoods(ensi) = likelihoods(ensi) + 0.5 * (H*xfc(:, l) - y(:, l))'*(R\(H*xfc(:, l) - y(:, l)));
                        end
                    end
                    
                    cumlikelihoods = [0; cumsum(likelihoods)];
                    %cumlikelihoods = cumlikelihoods - likelihoods(1);
                    
                    ensnchoosen = zeros(Nchoose, 1);
                    for j = 1:Nchoose
                        rv = cumlikelihoods(end)*rand();
                        zz = find(cumlikelihoods < rv);
                        
                        ensnchoosen(j) = zz(end);
                    end
                    
                    
                    ensM = xf(:, ensnchoosen, cea);
                    
                    analysis = mean(ensM, 2);
                    
                    %modelmodel.Y0 = mean(ensM(:, :, end), 2);
                    
                    % kill all humans
                    % TODO
                    
                    % calculate error
                    
                    erra = rms(analysis - xt(:, cea));
                    
                    err(i) = erra;
                    
                    errdis = err.^2;
                    if i <= spinup
                        errdis = nan;
                    else
                        errdis = sqrt(mean(errdis((spinup+1):i)));
                    end
                    
                    
                    fprintf(repmat('\b', 1, numel(dispstr)));
                    dispstr = sprintf('Err at it:%d, is: %.5f\n', i, errdis);
                    fprintf(dispstr);
                    
                end
                
                
                % propogate the truth and the model
                [~, Y] = solvernature(modeltruth.F, [0 dt], modeltruth.Y0);
                modeltruth.Y0 = Y(end, :).';
                
                xfprev = xf;
                
                
                % throw out
                if size(xf, 3) > L
                    xf = xf(:, :, 2:end);
                end
                
                ns = size(xf, 3) + 1;
                
                for ensi = 1:N
                    [~, Y] = ti(modelmodel.F, [0 dt], xf(:, ensi));
                    xf(:, ensi, ns) = Y(end, :).';
                end
                modelmodel.Y0 = mean(xf(:, :, end), 2);
            end
            
        end
    end
end
end


function [c, g] = cost(y, qens, L, ti, dt, F, Javp)



n = numel(y)/(L + 1);

y = reshape(y, n, L + 1);

e = zeros(L, 1);

d = zeros(n, L);

for l = 1:L
    [~, Y] = ti(F, [0 dt], y(:, l));
    mylp1 = Y(end, :).';
    
    d(:, l) = mylp1 + qens(:, l) - y(:, l + 1);
    
end

for l = 1:L
    e(l) = 0.5 * norm(d(:, l))^2;
end

c = sum(e);

if nargout > 1
    g = zeros(size(y));
    
    for l = 1:(L + 1)
        if l == 1
            %g(:, l) = d(:, l) + dt*Javp(0, y(:, l), d(:, l));
            g(:, l) = rk4adj(F, 0, y(:, l), Javp, dt, d(:, l));
        elseif l == (L + 1)
            g(:, l) = - d(:, l - 1);
        else
            %g(:, l) = - d(:, l - 1) + (d(:, l) + dt*Javp(0, y(:, l), d(:, l)));
            g(:, l) = rk4adj(F, 0, y(:, l), Javp, dt, d(:, l)) - d(:, l - 1);
        end
    end
    g = reshape(g, n*(L + 1), 1);
end

end


function flag = stop(x, y, H, R, L, sdallow)

n = numel(x)/(L + 1);

x = reshape(x, n, L + 1);

cs = 0;

for l = 1:(L + 1)
%     cs = cs + ((H*x(:, l) - y(:, l)).'*(R\(H*x(:, l) - y(:, l))))/(y(:, l).'*(R\y(:, l)));
cs = cs + ((H*x(:, l) - y(:, l)).'*(R\(H*x(:, l) - y(:, l))));
end

flag = sqrt(cs/(L + 1)) > sdallow;

end


function [t, y] = fe(f, t, y, steps)

tc = t(1);
dt = diff(t)/steps;

for step = 1:steps
    y = y + dt*f(tc, y);
    tc = tc + dt;
end

y = y.';

end


function [t, y] = rk4(f, t, y, dt)

tc = t(1);
k1 = dt*f(tc       , y       );
k2 = dt*f(tc + dt/2, y + k1/2);
k3 = dt*f(tc + dt/2, y + k2/2);
k4 = dt*f(tc + dt  , y + k3  );

y = y + k1/6 + k2/3 + k3/3 + k4/6;

y = y.';


end


function avp = rk4adj(f, t, y, Javp, dt, v)

tc = t(1);
k1 = dt*f(tc       , y       );
k2 = dt*f(tc + dt/2, y + k1/2);
k3 = dt*f(tc + dt/2, y + k2/2);
%k4 = dt*f(tc + dt  , y + k3  );

% dj1vp = dt*Javp(tc, y, v);
% dk1vp = dj1vp;
% 
% dj2vp = dt*Javp(tc, y + k1/2, v);
% dk2vp = dj2vp + 1/2*dt*Javp(tc, y, dj2vp);
% 
% dj3vp = Javp(tc, y + k2/2, v);
% jk1k2avp = Javp(tc, y + k1/2, dj3vp);
% j0k1k2avp = Javp(tc, y, jk1k2avp);
% dk3vp = dt*(dj3vp + 1/2*dt*jk1k2avp + 1/4*(dt^2)*j0k1k2avp);
% 
% dj4vp = Javp(tc, y + k3, v);
% jk2k3avp = Javp(tc, y + k2/2, dj4vp);
% jk1k2k3avp = Javp(tc, y + k1/2, jk2k3avp);
% j0k1k2k3avp = Javp(tc, y, jk1k2k3avp);
% dk4vp = dt*(dj4vp + 1/2*dt*jk2k3avp + 1/4*(dt^2)*jk1k2k3avp + 1/4*(dt^2)*j0k1k2k3avp);

%avp = v + dk1vp/6 + dk2vp/3 + dk3vp/3 + dk4vp/4;

lambda = v;

theta4 = dt*Javp(tc + dt  , y + k3  , (1/6*lambda)             );
theta3 = dt*Javp(tc + dt/2, y + k2/2, (1/3*lambda + 1*theta4));
theta2 = dt*Javp(tc + dt/2, y + k1/2, (1/3*lambda + 1/2*theta3));
theta1 = dt*Javp(tc       , y       , (1/6*lambda + 1/2*theta2));

avp = lambda + theta1 + theta2 + theta3 + theta4;


end
