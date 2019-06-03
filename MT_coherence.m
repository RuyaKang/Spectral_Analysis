% % Multitaper AR(4) and AR(2) coherence % %
AR4model = arima('Constant',0,'AR',{2.7607,-3.8106,2.6535,-0.9238},'Variance',1);
AR2model = arima('Constant', 0, 'AR', {0.75,-0.5}, 'Variance', 1);
% number of realisations
N = 1024;
% sampling interval
Delta = 1;
% bandwidth
NW = 8;

% frequency in [0, 0.5]
f = linspace(0,1/2,N/2+1);
% dpss up to order 2NW
dps = dpss(N, NW);

% obeservations in AR(4) process
X = simulate(AR4model,N);
% obeservations in AR(2) process
Y = simulate(AR2model,N);

% taper time series
hX = dps.*X;
hY = dps.*Y;

dseX = [];
dseY = [];
csdeXY = [];

for i = 1:2*NW
    dseX(i,:) = J(hX(:,i),Delta,f).*conj(J(hX(:,i),Delta,f));
    dseY(i,:) = J(hY(:,i),Delta,f).*conj(J(hY(:,i),Delta,f));
    csdeXY(i,:) = J(hX(:,i),Delta,f).*conj(J(hY(:,i),Delta,f));
end

% apply multitaper
dseXmt = transpose(1./(1:2*NW)).*cumsum(dseX);
dseYmt = transpose(1./(1:2*NW)).*cumsum(dseY);
csdeXYmt = transpose(1./(1:2*NW)).*cumsum(csdeXY);

% find coherence
figure()
gammaHat = [];
for i = 1:2*NW 
    gammaHat(i,:) = abs(csdeXYmt(i,:)).^2./(dseXmt(i,:).*dseYmt(i,:));
    subplot(NW,2,i)
    plot(f,gammaHat(i,:),'k');
    xlim([0 0.5]);
    ylim([0 1]);
    xlabel('f');
    ylabel(sprintf('K=%d',i)) 
end

% obeservations in combined AR(4) and AR(2) processes
W = X+Y;
Z = X-Y;

% (N*1) vector of sdf for frequency in [0, 0.5]
SX = ARsdf(4,[2.7607;-3.8106;2.6535;-0.9238],transpose(f));
SY = ARsdf(2,[0.75;-0.5],transpose(f));

% true coherence
gamma = abs(SX-SY).^2./(SX+SY).^2;

% taper time series
hW = dps.*W;
hZ = dps.*Z;

dseW = [];
dseZ = [];
csdeWZ = [];

for i = (1:2*NW)
    dseW(i,:) = J(hW(:,i),Delta,f).*conj(J(hW(:,i),Delta,f));
    dseZ(i,:) = J(hZ(:,i),Delta,f).*conj(J(hZ(:,i),Delta,f));
    csdeWZ(i,:) = J(hW(:,i),Delta,f).*conj(J(hZ(:,i),Delta,f));
end

% apply multitaper
dseWmt = transpose(1./(1:2*NW)).*cumsum(dseW);
dseZmt = transpose(1./(1:2*NW)).*cumsum(dseZ);
csdeWZmt = transpose(1./(1:2*NW)).*cumsum(csdeWZ);

% find coherence:,i
figure()
gammaHat = [];
for i = (1:2*NW)
    gammaHat(i,:) = abs(csdeWZmt(i,:)).^2./(dseWmt(i,:).*dseZmt(i,:));
    subplot(NW,2,i)
    plot(f,gammaHat(i,:),'k');
    hold on
    plot(f,gamma,'r');
    hold off
    xlim([0 0.5]);
    ylim([0 1]);
    xlabel('f');
    ylabel(sprintf('K=%d',i))
end
