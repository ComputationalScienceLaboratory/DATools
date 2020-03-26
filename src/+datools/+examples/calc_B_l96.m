m = otp.lorenz96.presets.Canonical;
[t, y] = ode45(m.Rhs.F, [0 10], m.Y0);

N = 500;

X = y(end, :).' + 1*randn(m.NumVars, N);

% prop

steps = 1000;

skip = 100;

B = zeros(m.NumVars);

for step = 1:steps
    step
    
    for i = 1:N
        [t, y] = datools.utils.rk4(m.Rhs.F, [0 0.05], X(:, i), 1);
        
        X(:, i) = y(end, :).';
    end
    
    
    if step > skip
        xfm = mean(X, 2);
        Af = (X - xfm)/sqrt(N - 1);
        
        B = B + (Af*Af.');
        imagesc(B/(step - skip)); colorbar; drawnow;
    end
end

B = B/(steps - skip);

save('l96B.mat', 'B');
