function LF = oracle(locfun, radbounds, nrs)

if nargin < 1
    locfun = @(~) 1;
end


LF = @(obj) locmatrix(obj, locfun, radbounds, nrs);

end


function rho = locmatrix(obj, locfun, radbounds, nrs)

global DAtruth;

R       = obj.CurrentObsErr;
ensN    = obj.EnsembleSize;
numVars = obj.Problem.NumVars;

Ys = obj.CurrentSyntheticObs;

%obsN = size(Ys,1);

H = obj.H;

% our current best guess is the forecast that we made before.
Xf =  obj.CurrentBestGuess;

Xfm = mean(Xf,2);
tmp = Xf-Xfm*ones(1,ensN);
Pf = 1/(ensN) *(tmp*tmp');



rPf = @(r) (locfun(numVars,r) .* Pf);
Hres     = @(A) A(H,:);
Hrest    = @(A) A(:,H);
Hresboth = @(A) A(H,H);


% S matrix from Law et al.
S = @(r) Hresboth(rPf(r))+R;

% Kalman gain matrix
K = @(r) (Hrest(rPf(r)))/S(r);

% Hey, analysis step of EnKF
Xa = @(r) Xf + K(r)*( Ys - Xf(H,:));

optFun = @(r) rms(mean(Xa(r),2) - DAtruth);

rmin = radbounds(1);
rmax = radbounds(end);

rmins = radbounds(1) * ones(nrs, 1);
rmaxs = radbounds(end) * ones(nrs, 1);


    

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
    
    rs = (rmins + rmaxs)/2;
    
    passes = 2;
    
    plot(1:numel(rs), rs, '-o');
    axis([0 numel(rs)+1 rmin-1 rmax+1]);
    title('Pass=1, ri=1');
    drawnow;
    
    rsperm = randperm(numel(rs));
    
    for pass = 1:passes
        for ii = rsperm
            replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
            Jnn = @(nr) optFun(replace(rs, ii, nr));
            allrs = rmin:0.25:rmax;
            
            [~, I] = min(arrayfun(Jnn, allrs));
            rs = replace(rs, ii, allrs(I));
            plot(1:numel(rs), rs, '-o');
            axis([0 numel(rs)+1 rmin-1 rmax+1]);
            title(sprintf('Pass=%d, ri=%d', pass, ii));
            drawnow;
        end
    end
    
    
    %% END



%fprintf('oracle r=%.5f|a=%.5f,aE=%.5f ||| r=5|a1.45,aE:%.5f\n',optr,opta,optFun(optra),optFun([5 1.45]));

rho = locfun(numVars, rs);

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