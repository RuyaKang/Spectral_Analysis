% % Spectral window % %
f = linspace(0,1/2,513);
t = linspace(1,32,32);
h = hann(32);
H = 10*log10(abs(sum(transpose(h).*exp(-j*2*pi*transpose(f).*t),2)).^2);
figure()
plot(f,H,'k')
xlim([0,1/2])
ylim([-40,20])
xlabel('f')
ylabel('H(f)')
