% % Fejer's kernel % %
f = linspace(0,1/2,513);
t = linspace(1,32,32);
F = 10*log10(1/32*abs(sum(exp(-j*2*pi*transpose(f).*t),2)).^2);
figure()
plot(f,F,'k')
xlim([0,1/2])
ylim([-40,20])
xlabel('f')
ylabel('F(f)')
