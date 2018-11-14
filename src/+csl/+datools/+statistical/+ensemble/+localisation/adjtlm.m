function LF = adjtlm(adjSolver, tlmSolver, ensembleMembers)

if nargin < 3
   ensembleMembers = inf; 
end

LF = @(obj) locmatrix(obj.Problem.NumVars, ensembleMembers, obj.Problem.F, adjSolver, tlmSolver, ...
    obj.PreviousTimeSpan, obj.PreviousAnalysis);

end


function rho = locmatrix(numVars, ensembleMembers, F, adjSolver, tlmSolver, previousTs, previousXa)
lambda = eye(numVars);

%xa = mean(previousXa,2);

rho = zeros(numVars);

ensN = size(previousXa,2);

if ensembleMembers > ensN
   ensembleMembers = ensN; 
end

ensNs = 1:ensembleMembers;

for j = ensNs
    
    xa = previousXa(:,j);
    
    [ ~, ~, adjSens ] = adjSolver(F,previousTs,xa,lambda);
    
    % here we have to restrict the set of vars whose TLM solution we will look
    % at for each variable.
    
    [ ~, ~, tlmSens ] = tlmSolver(F,previousTs,xa,adjSens);
    
    rho = rho + tlmSens;
    
end

%L = (1/ensN)*L;

rho = abs(rho);
%rho = abs(Pf);
% Let's normalise by the variance
mrho = diag(rho);
% so let's normalise by column.
%mL = max(L);
for k = 1:length(rho)
    rho(:,k) = rho(:,k)./mrho(k);
end

end