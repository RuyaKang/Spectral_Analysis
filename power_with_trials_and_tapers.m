alpha = 0.05;
Delta = 1;
NW = 5;

v = [0.5;0.5];
A = [0.3, 0.6; 0.6, 0.3];
B = [0.9, 0.8; 0.8, 0.9];
lambda0 = [0.5,0.5];

power = zeros(5,NW-1);
for i = 1:100
    sumSX = zeros(NW-1,51);
    sumSY = zeros(NW-1,51);
    sumSXY = zeros(NW-1,51);
    for n = 1:5
        [N,X,E,M,E_source,lambda,spectra,cross_spectra,coherences] = MvHawkesSimulation(v,A,B,lambda0);

        Xt = transpose(X(1,:));
        Yt = transpose(X(2,:));
        
        SX = [];
        SY = [];
        SXY = [];
        
        for k = 2:NW
            [SX(k-1,:),SY(k-1,:),SXY(k-1,:)] = sdf(Xt,Yt,Delta,k);
        end
        
        sumSX = sumSX+SX;
        sumSY = sumSY+SY;
        sumSXY = sumSXY+SXY;
        
        for j = 1:NW-1
            gammaHat = abs(sumSXY(j,:)/n).^2./((sumSX(j,:)/n).*(sumSY(j,:)/n));
            % rejection rate when goes through the stepdown procedure
            rate = mht(gammaHat,(k-1)*2*n,2,alpha);
            
            if rate==1
                power(n,j) = power(n,j)+1
            end
        end
        [i n]
    end
end

figure()
plot(1:5,power(:,1),'-o')
hold on
plot(1:5,power(:,2),'-o')
hold on
plot(1:5,power(:,3),'-o')
hold on
plot(1:5,power(:,4),'-o')
hold off
xlabel('Number of trials, k')
ylabel('Power (%)')
legend('NWt=2','NWt=3','NWt=4','NWt=5','Location','southeast')


