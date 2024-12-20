% this determines the distance between states in Lorenz96 used for
% determining the distance matrix
%
% this is exactly in accordance to OTP 
% reference (?)
function d = distanceFunc(~, y, i, j)
% i, j are the states concerned
N = numel(y);
d = min([abs(i-j); abs(N+i-j); abs(N+j-i)]);
end