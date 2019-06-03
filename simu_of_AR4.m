% % AR(4) realisation % %
AR4model = arima('Constant',0,'AR',{2.7607,-3.8106,2.6535,-0.9238},'Variance',1);
X = simulate(AR4model,1024);

figure()
plot(1:1024,X,'k')
xlim([1 1024])
xlabel('t')
ylabel('Xt')
