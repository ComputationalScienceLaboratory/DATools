function [alpha, steps] = armijowolfe(ffun, gfun, x, p, g)

%% Constant setup
alphahat = 1;
c1 = 1e-4;
c2 = 0.9;
alphamin = 1e-3;
maxbacklineits = 10;
%maxbacklineits = 50;

%% Initial setup
alpha    = alphahat;
Phi0     = ffun(x);
Phihat   = ffun(x + alphahat*p);
dPhi0    = p.' * g;
Phialpha = Phihat;
wolfeg   = gfun(x + alpha*p);

steps = 0;

% Check the Armijo, and the strong Wolfe conditions

while (Phialpha > Phi0 + c1*alpha*dPhi0 || ...
        abs(p.' * wolfeg) > c2*abs(dPhi0)) && ...
        alpha >= alphamin && ...
        steps < maxbacklineits
    
    
    if steps > 0
        %% Cubic interpolation
        d = Phi0;
        c = dPhi0;
        a = -(-d*alpha^2 + Phihat*alpha^2 - c*(alpha^2)*alphahat + ...
            d*alphahat^2 - Phialpha*alphahat^2 + c*alpha*(alphahat^2))/ ...
            ((alpha^2)*(alpha - alphahat)*(alphahat^2));
        b = -(d*(alpha^3) - Phihat*(alpha^3) + c*(alpha^3)*alphahat - ...
            d*(alphahat^3) + Phialpha*(alphahat^3) - c*alpha*(alphahat^3))/ ...
            ((alpha^2)*(alpha - alphahat)*(alphahat^2));
        
        disc = sqrt(b^2 - 3*a*c);
        
        alphaplus  = (-b + disc)/(3*a);
        alphaminus = (-b - disc)/(3*a);
        
        Phihalphaplus = ffun(x + alphaplus*p);
        Phihalphaminus = ffun(x + alphaminus*p);

        if Phihalphaplus < Phihalphaminus
            alpha = alphaplus;
            Phialpha = Phihalphaplus;
        else
            alpha = alphaminus;
            Phialpha = Phihalphaminus;
        end  
        
        
    else
        %% Quadratic interpolation
        alpha = -(dPhi0*alphahat)/(2*(Phihat - Phi0 - dPhi0*alphahat));
        Phialpha = ffun(x + alpha*p);
    end
    

    wolfeg   = gfun(x + alpha*p);
    steps = steps + 1;
end

%% Heuristic
alpha = max(alpha, alphamin);

end
