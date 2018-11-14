function lambda = rk4adj(f, t, y, Javp, lambda)

for ti = (numel(t) - 1):-1:1

    yc = y(ti, :).';
    
    tc = t(ti);
    
    dt = t(ti + 1) - tc;
    
    k1 = dt*f(tc       , yc       );
    k2 = dt*f(tc + dt/2, yc + k1/2);
    k3 = dt*f(tc + dt/2, yc + k2/2);
    
    theta4 = dt*Javp(tc + dt  , yc + k3  , (1/6*lambda)             );
    theta3 = dt*Javp(tc + dt/2, yc + k2/2, (1/3*lambda + 1*theta4));
    theta2 = dt*Javp(tc + dt/2, yc + k1/2, (1/3*lambda + 1/2*theta3));
    theta1 = dt*Javp(tc       , yc       , (1/6*lambda + 1/2*theta2));
    
    lambda = lambda + theta1 + theta2 + theta3 + theta4;

end

end
