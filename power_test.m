alpha = 0.05;
Delta = 1;
NW = 4;

v = [0.5;0.5];
B = [0.9, 0.8; 0.8, 0.9];
lambda0 = [0.5,0.5];

power = zeros(2,5);

for i = 1:5
    A = [0.3,0.4+0.1*i;0.4+0.1*i,0.3];
    for j = 1:100
        [N,X,E,M,E_source,lambda,spectra,cross_spectra,coherences] = MvHawkesSimulation(v,A,B,lambda0);
        Xt = transpose(X(1,:));
        Yt = transpose(X(2,:));
        
        gammaHat = coh(Xt,Yt,NW,Delta);
        rate = mht(gammaHat,(NW-1)*2,2,alpha);
        
        if rate==1
            power(1,i) = power(1,i)+1
        end
        
        C = fzero(@(x) estcoh(x,(NW-1)*2,2,(1-alpha)^(1/length(Xt)),0), [0.000001, 0.999999]);
        
        if min(gammaHat)>C
            power(2,i) = power(2,i)+1;
        end
        [i j]
    end
end

figure()
plot(0.4+0.1*[1:5],power(1,:),'-o')
hold on
plot(0.4+0.1*[1:5],power(2,:),'-o')
hold off
xlabel('Mutual excitation parameter')
ylabel('Power (%)')
legend('Stepdown','Standard Hypithesis test','Location','southeast')

