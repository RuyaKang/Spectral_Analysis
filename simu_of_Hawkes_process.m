% % Hawles process realisation % %
[test_counts,test_times,test_lambda,test_theory,test_spectrum] = st_hawkes_1d_simulator(0.5,0.8,1.2,0.5);
subplot(2,1,1)
plot(test_times,cumsum(histcounts(test_times,0:max(test_times)/length(test_times):max(test_times))),'k')
xlim([0 max(test_times)]) 
xlabel('t')
ylabel('Count')
subplot(2,1,2)
plot(test_times,test_lambda,'k')
xlim([0 max(test_times)])
xlabel('t')
ylabel('Intensity')
