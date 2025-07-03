% this is the jacobian of the rhs of the lorenz96 equations
%
function J = jacobian(~, y, N)
% J = zeros(N,N);
mainDiag = -1*ones(N,1);
dPlusOne = zeros(N-1,1);
dMinusOne = zeros(N-1,1);
dminusTwo = zeros(N-2,1);

for i = 1:N-1
    dPlusOne(i) = y(mod(i - 1 - 1, N) + 1);
end

for i = 2:N
    dMinusOne(i) = - y(mod(i - 2 - 1, N) + 1);
end

for i = 3:N
    dminusTwo(i) = - y(mod(i - 1 - 1, N) + 1);
end

J = diag(mainDiag) + diag(dPlusOne) + diag(dMinusOne) + diag(dminusTwo);

end