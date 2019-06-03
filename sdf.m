function [SXmt,SYmt,SXYmt] = sdf(X,Y,Delta,NW)

% multitaper spectrum
N = length(X);

f = linspace(0,1/(2*Delta),N+1);
% dps sequences up to order 2NW
dps = dpss(N,NW);

% taper time series
hX = dps.*X;
hY = dps.*Y;

dseX = [];
dseY = [];
csdeXY = [];

for k = (1:(NW-1)*2)
    dseX(k,:) = J(hX(:,k),Delta,f).*conj(J(hX(:,k),Delta,f));
    dseY(k,:) = J(hY(:,k),Delta,f).*conj(J(hY(:,k),Delta,f));
    csdeXY(k,:) = J(hX(:,k),Delta,f).*conj(J(hY(:,k),Delta,f));
end

SXmt = mean(dseX,1);
SYmt = mean(dseY,1);
SXYmt = mean(csdeXY,1);