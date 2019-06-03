function J = J(X, Delta, f)
% compute J for a range of frequency
%     input: 
%        X = column vector of observations
%        Delta = sampling interval
%        f = row vector of frequencies
%     output:
%        the defined function J

% number of observations
N = length(X);
% J(f) for the data taper
J = sqrt(Delta) * exp(-j*2*pi*transpose(f).*(1:N)*Delta)*X;
end