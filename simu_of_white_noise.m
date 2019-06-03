% % White noise realisation with mean 0 and variance 1 % %
X = [];
for i = 1:1024
    X(i) = normrnd(0,1);
end

figure()
plot(1:1024,X,'k')
xlim([1 1024])
xlabel('t')
ylabel('Xt')