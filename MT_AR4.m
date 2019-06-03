
% % Multitaper AR(4) sdf % %
AR4model = arima('Constant',0,'AR',{2.7607,-3.8106,2.6535,-0.9238},'Variance',1);
% number of realisations
N = 1024;
% sampling interval
Delta = 1;
% bandwidth
NW = 4;

% frequency in [0, 0.5]
f = linspace(0,1/2,N/2+1);
% obeservations in AR(4) process
X = simulate(AR4model, N);
% sdf for frequency in [0, 0.5] in dB
S = 10*log10(ARsdf(4,[2.7607;-3.8106;2.6535;-0.9238],transpose(f)));

% original AR(4) process and its sdf
subplot(2,1,1)
plot(1:N,X,'k');
xlim([0 1024]);
xlabel('t');
ylabel('Xt')

subplot(2,1,2)
plot(f,S,'k');
xlim([0 0.5]);
xlabel('f');
ylabel('S(f) in dB');

% dpss up to order 2NW
dps= dpss(N,NW);
% taper time series
hX = dps.*X;

% taper and tapered AR(4) process for k = 0 to 3
figure()
for i = 1:4  
    subplot(4,4,4*i-3)
    plot(1:N,dps(:,i),'k');
    xlim([0,1024]);
    ylim([-0.1,0.1]);
    xlabel('t');
    ylabel(sprintf('k=%d',i-1));
    
    subplot(4,4,4*i-2)
    plot((1:N),hX(:,i),'k');
    xlim([0 1024]);
    xlabel('t');

end

% taper and tapered AR(4) process for k = 4 to 7
for i = 5:8
    subplot(4,4,4*(i-4)-1)
    plot(1:N,dps(:,i),'k');
    xlim([0 1024]);
    xlabel('t');
    ylabel(sprintf('k=%d',i-1));
    
    subplot(4,4,4*(i-4))
    plot((1:N),hX(:,i),'k');
    xlim([0 1024]);
    xlabel('t');
end

% decrease the range of frequency to get a better plot for spectral windows
f1 =  linspace(0, 0.01, N);
% spectral window and direct spectral estimates in dB
H = [];
dseX = [];
for i = 1:8
    H(i,:) = J(dps(:,i),Delta,f1).*conj(J(dps(:,i),Delta,f1));
    dseX(i,:) = J(hX(:,i),Delta,f).*conj(J(hX(:,i),Delta,f));
end

% apply multitaper
Hmt = transpose(1./(1:2*NW)).*cumsum(H);
dsemt = transpose(1./(1:2*NW)).*cumsum(dseX);

% mean of N spectral window and periodogram for N = 1 to 4 in dB
figure()
for i = 1:4
    subplot(4,4,4*i-3)
    plot(f1,10*log10(Hmt(i,:)),'k');
    xlim([0 0.01]);
    ylim([-80 40]);
    xlabel('f');
    ylabel('H(f)');
    title(sprintf('K=%d',i))
    
    subplot(4,4,4*i-2)
    plot(f,10*log10(dsemt(i,:)),'k');
    hold on 
    plot(f,S,'r');
    hold off
    xlim([0 0.5]);
    ylim([-40 60]);
    xlabel('f');
    ylabel('S(f)');
    title(sprintf('K=%d',i))

end

% mean of N spectral window and periodogram for N = 5 to 8 in dB
for i = 5:8
    subplot(4,4,4*(i-4)-1)
    plot(f1,10*log10(Hmt(i,:)),'k');
    xlim([0 0.01]);
    ylim([-80 40]);
    xlabel('f');
    ylabel('H(f)');
    title(sprintf('K=%d',i))
    
    subplot(4,4,4*(i-4))
    plot(f,10*log10(dsemt(i,:)),'k');
    hold on 
    plot(f,S,'r');
    hold off
    xlim([0 0.5]);
    ylim([-40 60]);
    xlabel('f');
    ylabel('S(f)');
    title(sprintf('K=%d',i))
end
