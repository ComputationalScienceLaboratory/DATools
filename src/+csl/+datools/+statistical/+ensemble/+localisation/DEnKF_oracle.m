function LF = DEnKF_oracle(locfun, radbounds, nrs, distfn, optr, m)

if nargin < 1
    locfun = @(~) 1;
end

if nargin < 6
    m = @(ri, rj) (ri + rj)/2;
end


LF = @(obj, H) locmatrix(obj, locfun, radbounds, nrs, H, distfn, m, optr);

end


function rho = locmatrix(obj, locfun, radbounds, nrs, H, distfn, m, optr)

R       = obj.CurrentObsErr;
ensN    = obj.EnsembleSize;
numVars = obj.Problem.NumVars;

Y = obj.CurrentObs;

%obsN = size(Ys,1);

% our current best guess is the forecast that we made before.
Xf =  obj.CurrentBestGuess;

Xfm = mean(Xf, 2);

d = Y - Xfm(H, :);

Af = Xf-Xfm*ones(1,ensN);
PfH = 1/(ensN) *(Af*Af(H,:).');


%rs = fminunc(@(r) cost(r, Xfm, PfH, d, H, R, numVars, locfun, distfn, radbounds, m), 0.5, psoptions);

beta = optr*ones(nrs, 1);


lb = radbounds(1)*ones(nrs, 1);
ub = radbounds(2)*ones(nrs, 1);

global DAtruth;
Jn = @(nb) cost(nb, Xfm, PfH, d, H, R, numVars, locfun, distfn, radbounds, m, DAtruth);

if nrs > 4
    psoptions = optimoptions('fminunc', 'Display', 'off', ...
        'UseParallel', true, 'FunctionTolerance', eps, 'OptimalityTolerance', eps);
    %particleswarm

else
    psoptions = optimoptions('fminunc', 'Display', 'off', ...
        'UseParallel', false, 'FunctionTolerance', eps, 'OptimalityTolerance', eps);
end

% fmuopts = optimoptions('fminunc', 'Display', 'off', ...
%  'UseParallel', false, 'Algorithm', 'quasi-newton');
rsperm = 1:nrs;
for ii = rsperm
    replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
    Jnn = @(nb) Jn(replace(beta, ii, nb));
    %fminunc(@(r) , 0.5, psoptions);
    zz = fminunc(Jnn, beta(ii), psoptions);
    
    %[~, I] = min(arrayfun(Jnn, allrs));
    beta = replace(beta, ii, zz);
end

%beta = patternsearch(Jn, beta, [],[],[],[],lb, ub,[], psoptions);

%beta = fminunc(Jn, beta, [],[],[],[],lb, ub,[], psoptions);
%beta = fminunc(Jn, beta, psoptions);

rs = beta(:);


% ndar = 50;
% radii = linspace(rmin,rmax,ndar);
% infs  = linspace(amin,amax,ndar);
% Jevals = zeros(numel(radii),numel(infs));
% for i = 1:numel(radii)
%     for j = 1:numel(infs)
%         Jevals(i,j) = optFun([radii(i),infs(j)]);
%     end
% end
% surf(infs,radii,Jevals);
% drawnow;


    %psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 1.5);
    %optr = patternsearch(optFun,rmin,[],[],[],[],rmin,rmax,[],psoptions);
    
    %% ALTERNATIVE
    
%     rs = (rmins + rmaxs)/2;
%     
%     passes = 1;
%     
%     plot(1:numel(rs), rs, '-o');
%     axis([0 numel(rs)+1 rmin-1 rmax+1]);
%     title('Pass=1, ri=1');
%     drawnow;
%     
%     rsperm = randperm(numel(rs));
%     
%     for pass = 1:passes
%         for ii = rsperm
%             replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
%             Jnn = @(nr) optFun(replace(rs, ii, nr));
%             allrs = rmin:0.25:rmax;
%             
%             [~, I] = min(arrayfun(Jnn, allrs));
%             rs = replace(rs, ii, allrs(I));
%             plot(1:numel(rs), rs, '-o');
%             axis([0 numel(rs)+1 rmin-1 rmax+1]);
%             title(sprintf('Pass=%d, ri=%d', pass, ii));
%             drawnow;
%         end
%     end
    
    
    %% END



%fprintf('oracle r=%.5f|a=%.5f,aE=%.5f ||| r=5|a1.45,aE:%.5f\n',optr,opta,optFun(optra),optFun([5 1.45]));

rho = locfun(numVars, rs, distfn, [], [], [], [], H);

end


function c = cost(r, Xfm, PfH, d, H, R, numVars, locfun, distfn, rb, m, DAtruth)


warning off;

r = r(:);

rPfH = locfun(numVars, r, distfn, [], [], m, [], H) .* PfH;

S = rPfH(H, :) + R;

K = rPfH/S;

Xam = Xfm + K*d;

c = norm(Xam - DAtruth);

if any(r > rb(2))
    c = inf;
end
if any(r < rb(1))
    c = inf;
end

end

% function rho = locfun(n,r)
% rho = gasparicohnmat(n,r);
% end
% 
% function L = gaussmat(n,r)
% L = zeros(n, n);
% for i = 1:n
%     for j = 1:n
%         d = min([abs(i-j),abs(n+i-j),abs(n+j-i)]);
%         L(i,j) = exp(-(d^2)/(2*r^2));
%     end
% end
% end
% 
% function L = gasparicohnmat(n,r)
% L = zeros(n, n);
% for i = 1:n
%     for j = 1:n
%         d = min([abs(i-j),abs(n+i-j),abs(n+j-i)]);
%         L(i,j) = gc(d/r,1);
%     end
% end
% end
% 
% function g = gc(k,theta)
% if (k <= theta)
%     g = 1 - (5/3)*k^2 + (5/8) * k^3 + (1/2)*k^4 - (1/4)*k^5;
% elseif (k <= theta)
%     g = 4 - 5*k+ (5/3)*k^2+(5/8)*k^3-(1/2)*k^4+(1/12)*k^5-(2/(3*k));
% else
%     g = 0;
% end
% end