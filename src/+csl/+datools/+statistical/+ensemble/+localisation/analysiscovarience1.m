function LF = analysiscovarience1
LF = @(obj) analysiscovariencemat(obj);
end


function rho = analysiscovariencemat(obj)
R = obj.CurrentObsErr;

numVars = obj.Problem.NumVars;


H = obj.H;


ensN = obj.EnsembleSize;
Xf =  obj.CurrentBestGuess;

numrounds = 1;

rho = ones(numVars);

for i = 1:numrounds
    
    % We then calculate the mean of all the ensemble members
    Xfm = mean(Xf,2);
    
    % And the covarience
    tmp = Xf-Xfm*ones(1,ensN);
    Pf = 1/(ensN) *(tmp*tmp');
    
    % Shur Product
    Pf = Pf .* rho;
    
    % S matrix from Law et al.
    S = Pf(H,H)+R;
    
    % Kalman gain matrix
    K = (Pf(:,H))/S;
    
    Ys = obj.CurrentSyntheticObs;
    
    % Hey, analysis step of EnKF
    %Xf = Xf - K*( Xf(H,:) - Ys);
    
    %Xfm = mean(Xf,2);
    %tmp = Xf-Xfm*ones(1,ensN);
    %rho = 1/(ensN) *(tmp*tmp');
    
    rho = Pf;
    
    rho = abs(rho);
    %rho = abs(Pf);
    % Let's normalise by the variance
    %mrho = diag(rho);
    % so let's normalise by column.
    mrho = max(rho);
    for k = 1:length(rho)
        rho(:,k) = rho(:,k)./mrho(k);
    end
end

end
