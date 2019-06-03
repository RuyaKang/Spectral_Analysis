Delta = 0.14;
NW = 4;
alpha = 0.05;

for sfi = 1:10
    
    rate = zeros(185,185);
    
    for m = 1:185
        for n = m+1:185
            sumSX = 0;
            sumSY = 0;
            sumSXY = 0;
            for i = 1:6
                sfc1 = get_sf(cells, m,i+6*(sfi-1));
                proc1_freq1 = histcounts(sfc1,1:Delta:8);
                sf1 = transpose(proc1_freq1-mean(proc1_freq1));

                sfc2 = get_sf(cells,n,i+6*(sfi-1));
                proc1_freq2 = histcounts(sfc2,1:Delta:8);
                sf2 = transpose(proc1_freq2-mean(proc1_freq2));

                [SX,SY,SXY] = sdf(sf1,sf2,Delta,NW);
                sumSX = sumSX+SX;
                sumSY = sumSY+SY;
                sumSXY = sumSXY+SXY;
            end

            gammaHat = abs(sumSXY).^2./(sumSX.*sumSY);
            % rejection rate when goes through the stepdown procedure
            rate(m,n) = mht(gammaHat,(NW-1)*2*4,2,alpha);
            [m n]
        end
    end

    gamma_mat = rate;

    for i = 1:185
        for j = 1:i
            if i==j
                gamma_mat(i,j)=1;
            else
                gamma_mat(i,j) = gamma_mat(j,i);
            end
        end
    end
    
    figure()
    heatmap(gamma_mat)

end