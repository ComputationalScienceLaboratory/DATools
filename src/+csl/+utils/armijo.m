function alpha = Armijo(fhandle, x, p, g)

c1 = 1e-4;

alpha = 1;

Phi0 = fhandle(x);
Phi1 = fhandle(x + p);
dPhi0 = p.' * g;
Phialpha = Phi1;

alphamin = 1e-3;

steps = 0;

while Phialpha > Phi0 + c1*alpha*dPhi0 && alpha >= alphamin
    
    if steps > 1 && alpha < 1
        d = Phi0;
        c = dPhi0;
        a = -(d - Phialpha + c*alpha - ...
            c*alpha^2 - d*alpha^2 + Phi1*alpha^2)/...
            (-alpha^2 + alpha^3);
        b = -(-d + Phialpha - c*alpha + ...
            c*alpha^3 + d*alpha^3 - Phi1*alpha^3)/...
            (-alpha^2 + alpha^3);
        
        disc = sqrt(b^2 - 3*a*c);
        
        alphaplus  = (-b + disc)/(3*a);
        alphaminus = (-b - disc)/(3*a);

        if fhandle(x + alphaplus*p) < fhandle(x + alphaminus*p)
            alpha = alphaplus;
        else
            alpha = alphaminus;
        end
        
    else
        alpha = -(dPhi0)/(2*(Phi1 - Phi0 - dPhi0));
    end
    
    Phialpha = fhandle(x + alpha*p);   
    
    steps = steps + 1;
end

alpha = max(alpha, alphamin);

end
