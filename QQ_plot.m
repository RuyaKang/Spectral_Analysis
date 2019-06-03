% number of trials to generate the QQ-plot
trial = 100;

% basic settings
AR4model = arima('Constant', 0, 'AR', {2.7607, -3.8106, 2.6535, -0.9238}, 'Variance', 1);
AR2model = arima('Constant', 0, 'AR', {0.75, -0.5}, 'Variance', 1);
N = 128;
% sampling interval
delta_t = 1;
% bandwidth
w = 4/N;
% number of tapers
K = 4;
% the choice of frequency
n = 3;

% (N/2*1) vector frequency in [0, 0.5]
f = linspace(0,1/2,N+1);
% (2NW*N) matrix of dps sequences up to order 2NW
dps = dpss(N, N*w);

% (N*1) vector of sdf for frequency in [0, 0.5]
S_X = ARsdf(4, [2.7607; -3.8106; 2.6535; -0.9238], transpose(f(2:end)));
S_Y = ARsdf(2, [0.75; -0.5], transpose(f(2:end)));

% true coherence
gamma = abs(S_X-S_Y).^2./(S_X+S_Y).^2;

gamma1_f_WN = [];
gamma1_f = [];
gamma1_f_WZ = [];

for i = (1:trial)
    % (1*N) vector of obeservations in white noise processes with burn-in's
    WN1 = transpose(normrnd(0, 1, [1, N]));
    WN2 = transpose(normrnd(0, 1, [1, N]));
    % (1*N) vector of obeservations in AR(4) and AR(2) processes with burn-in's
    X = simulate(AR4model, N);
    Y = simulate(AR2model, N);
    % (1*N) vector of obeservations in AR(4) and AR(2) processes with burn-in's
    W = X + Y;
    Z = X - Y;

    % taper time series
    hWN1 = dps.*WN1;
    hWN2 = dps.*WN2;
    % taper independent time series
    hX = dps.*X;
    hY = dps.*Y;
    % taper time series
    hW = dps.*W;
    hZ = dps.*Z;
    
    dse1 = []; dse2 = []; csde12 = [];
    dse_X = []; dse_Y = []; csde_XY = [];
    dse_W = []; dse_Z = []; csde_WZ = [];
    
    for j = (1:2*N*w)
        dse1(j, :) = J(hWN1(:,j), delta_t, f).*conj(J(hWN1(:,j), delta_t, f));
        dse2(j, :) = J(hWN2(:,j), delta_t, f).*conj(J(hWN2(:,j), delta_t, f));
        csde12(j, :) = J(hWN1(:,j), delta_t, f).*conj(J(hWN2(:,j), delta_t, f));
        
        dse_X(j, :) = J(hX(:,j), delta_t, f).*conj(J(hX(:,j), delta_t, f));
        dse_Y(j, :) = J(hY(:,j), delta_t, f).*conj(J(hY(:,j), delta_t, f));
        csde_XY(j, :) = J(hX(:,j), delta_t, f).*conj(J(hY(:,j), delta_t, f));
        
        dse_W(j, :) = J(hW(:,j), delta_t, f).*conj(J(hW(:,j), delta_t, f));
        dse_Z(j, :) = J(hZ(:,j), delta_t, f).*conj(J(hZ(:,j), delta_t, f));
        csde_WZ(j, :) = J(hW(:,j), delta_t, f).*conj(J(hZ(:,j), delta_t, f));
    end
    
    % apply multitaper
    dse1_mt = transpose(1./(1:2*N*w)).*cumsum(dse1);
    dse2_mt = transpose(1./(1:2*N*w)).*cumsum(dse2);
    csde12_mt = transpose(1./(1:2*N*w)).*cumsum(csde12);
    % apply multitaper to AR processes
    dse_X_mt = transpose(1./(1:2*N*w)).*cumsum(dse_X);
    dse_Y_mt = transpose(1./(1:2*N*w)).*cumsum(dse_Y);
    csde_XY_mt = transpose(1./(1:2*N*w)).*cumsum(csde_XY);
    % apply multitaper to dependent processes
    dse_W_mt = transpose(1./(1:2*N*w)).*cumsum(dse_W);
    dse_Z_mt = transpose(1./(1:2*N*w)).*cumsum(dse_Z);
    csde_WZ_mt = transpose(1./(1:2*N*w)).*cumsum(csde_WZ);

    % find coherence
    gamma_hat_WN = [];
    gamma_hat = [];
    gamma_hat_WZ = [];
    
    for j = (1:2*N*w)
        gamma_hat_WN(j, :) = abs(csde12_mt(j, :)).^2./(dse1_mt(j, :).*dse2_mt(j, :));
        gamma_hat(j, :) = abs(csde_XY_mt(j, :)).^2./(dse_X_mt(j, :).*dse_Y_mt(j, :));
        gamma_hat_WZ(j, :) = abs(csde_WZ_mt(j, :)).^2./(dse_W_mt(j, :).*dse_Z_mt(j, :));
    end
    
    gamma1_f_WN(i) = gamma_hat_WN(K, n);
    gamma1_f(i) = gamma_hat(K, n);
    gamma1_f_WZ(i) = gamma_hat_WZ(K, n);
end

% number of trials to generate the QQ-plot
trial = 100;

% basic settings
AR4model = arima('Constant', 0, 'AR', {2.7607, -3.8106, 2.6535, -0.9238}, 'Variance', 1);
AR2model = arima('Constant', 0, 'AR', {0.75, -0.5}, 'Variance', 1);
N = 128;
% sampling interval
delta_t = 1;
% bandwidth
w = 4/N;
% number of tapers
K = 4;
% the choice of frequency
n = 50;

% (N/2*1) vector frequency in [0, 0.5]
f = linspace(0,1/2,N+1);
% (2NW*N) matrix of dps sequences up to order 2NW
dps = dpss(N, N*w);

% (N*1) vector of sdf for frequency in [0, 0.5]
S_X = ARsdf(4, [2.7607; -3.8106; 2.6535; -0.9238], transpose(f(2:end)));
S_Y = ARsdf(2, [0.75; -0.5], transpose(f(2:end)));

% true coherence
gamma = abs(S_X-S_Y).^2./(S_X+S_Y).^2;

gamma2_f_WN = [];
gamma2_f = [];
gamma2_f_WZ = [];

for i = (1:trial)
    % (1*N) vector of obeservations in white noise processes with burn-in's
    WN1 = transpose(normrnd(0, 1, [1, N]));
    WN2 = transpose(normrnd(0, 1, [1, N]));
    % (1*N) vector of obeservations in AR(4) and AR(2) processes with burn-in's
    X = simulate(AR4model, N);
    Y = simulate(AR2model, N);
    % (1*N) vector of obeservations in AR(4) and AR(2) processes with burn-in's
    W = X + Y;
    Z = X - Y;

    % taper time series
    hWN1 = dps.*WN1;
    hWN2 = dps.*WN2;
    % taper independent time series
    hX = dps.*X;
    hY = dps.*Y;
    % taper time series
    hW = dps.*W;
    hZ = dps.*Z;
    
    dse1 = []; dse2 = []; csde12 = [];
    dse_X = []; dse_Y = []; csde_XY = [];
    dse_W = []; dse_Z = []; csde_WZ = [];
    
    for j = (1:2*N*w)
        dse1(j, :) = J(hWN1(:,j), delta_t, f).*conj(J(hWN1(:,j), delta_t, f));
        dse2(j, :) = J(hWN2(:,j), delta_t, f).*conj(J(hWN2(:,j), delta_t, f));
        csde12(j, :) = J(hWN1(:,j), delta_t, f).*conj(J(hWN2(:,j), delta_t, f));
        
        dse_X(j, :) = J(hX(:,j), delta_t, f).*conj(J(hX(:,j), delta_t, f));
        dse_Y(j, :) = J(hY(:,j), delta_t, f).*conj(J(hY(:,j), delta_t, f));
        csde_XY(j, :) = J(hX(:,j), delta_t, f).*conj(J(hY(:,j), delta_t, f));
        
        dse_W(j, :) = J(hW(:,j), delta_t, f).*conj(J(hW(:,j), delta_t, f));
        dse_Z(j, :) = J(hZ(:,j), delta_t, f).*conj(J(hZ(:,j), delta_t, f));
        csde_WZ(j, :) = J(hW(:,j), delta_t, f).*conj(J(hZ(:,j), delta_t, f));
    end
    
    % apply multitaper
    dse1_mt = transpose(1./(1:2*N*w)).*cumsum(dse1);
    dse2_mt = transpose(1./(1:2*N*w)).*cumsum(dse2);
    csde12_mt = transpose(1./(1:2*N*w)).*cumsum(csde12);
    % apply multitaper to AR processes
    dse_X_mt = transpose(1./(1:2*N*w)).*cumsum(dse_X);
    dse_Y_mt = transpose(1./(1:2*N*w)).*cumsum(dse_Y);
    csde_XY_mt = transpose(1./(1:2*N*w)).*cumsum(csde_XY);
    % apply multitaper to dependent processes
    dse_W_mt = transpose(1./(1:2*N*w)).*cumsum(dse_W);
    dse_Z_mt = transpose(1./(1:2*N*w)).*cumsum(dse_Z);
    csde_WZ_mt = transpose(1./(1:2*N*w)).*cumsum(csde_WZ);

    % find coherence
    gamma_hat_WN = [];
    gamma_hat = [];
    gamma_hat_WZ = [];
    
    for j = (1:2*N*w)
        gamma_hat_WN(j, :) = abs(csde12_mt(j, :)).^2./(dse1_mt(j, :).*dse2_mt(j, :));
        gamma_hat(j, :) = abs(csde_XY_mt(j, :)).^2./(dse_X_mt(j, :).*dse_Y_mt(j, :));
        gamma_hat_WZ(j, :) = abs(csde_WZ_mt(j, :)).^2./(dse_W_mt(j, :).*dse_Z_mt(j, :));
    end
    
    gamma2_f_WN(i) = gamma_hat_WN(K, n);
    gamma2_f(i) = gamma_hat(K, n);
    gamma2_f_WZ(i) = gamma_hat_WZ(K, n);
end

% number of trials to generate the QQ-plot
trial = 100;

% basic settings
AR4model = arima('Constant', 0, 'AR', {2.7607, -3.8106, 2.6535, -0.9238}, 'Variance', 1);
AR2model = arima('Constant', 0, 'AR', {0.75, -0.5}, 'Variance', 1);
N = 1024;
% sampling interval
delta_t = 1;
% bandwidth
w = 4/N;
% number of tapers
K = 4;
% the choice of frequency
n = 24;

% (N/2*1) vector frequency in [0, 0.5]
f = linspace(0,1/2,N+1);
% (2NW*N) matrix of dps sequences up to order 2NW
dps = dpss(N, N*w);

% (N*1) vector of sdf for frequency in [0, 0.5]
S_X = ARsdf(4, [2.7607; -3.8106; 2.6535; -0.9238], transpose(f(2:end)));
S_Y = ARsdf(2, [0.75; -0.5], transpose(f(2:end)));

% true coherence
gamma = abs(S_X-S_Y).^2./(S_X+S_Y).^2;

gamma3_f_WN = [];
gamma3_f = [];
gamma3_f_WZ = [];

for i = (1:trial)
    % (1*N) vector of obeservations in white noise processes with burn-in's
    WN1 = transpose(normrnd(0, 1, [1, N]));
    WN2 = transpose(normrnd(0, 1, [1, N]));
    % (1*N) vector of obeservations in AR(4) and AR(2) processes with burn-in's
    X = simulate(AR4model, N);
    Y = simulate(AR2model, N);
    % (1*N) vector of obeservations in AR(4) and AR(2) processes with burn-in's
    W = X + Y;
    Z = X - Y;

    % taper time series
    hWN1 = dps.*WN1;
    hWN2 = dps.*WN2;
    % taper independent time series
    hX = dps.*X;
    hY = dps.*Y;
    % taper time series
    hW = dps.*W;
    hZ = dps.*Z;
    
    dse1 = []; dse2 = []; csde12 = [];
    dse_X = []; dse_Y = []; csde_XY = [];
    dse_W = []; dse_Z = []; csde_WZ = [];
    
    for j = (1:2*N*w)
        dse1(j, :) = J(hWN1(:,j), delta_t, f).*conj(J(hWN1(:,j), delta_t, f));
        dse2(j, :) = J(hWN2(:,j), delta_t, f).*conj(J(hWN2(:,j), delta_t, f));
        csde12(j, :) = J(hWN1(:,j), delta_t, f).*conj(J(hWN2(:,j), delta_t, f));
        
        dse_X(j, :) = J(hX(:,j), delta_t, f).*conj(J(hX(:,j), delta_t, f));
        dse_Y(j, :) = J(hY(:,j), delta_t, f).*conj(J(hY(:,j), delta_t, f));
        csde_XY(j, :) = J(hX(:,j), delta_t, f).*conj(J(hY(:,j), delta_t, f));
        
        dse_W(j, :) = J(hW(:,j), delta_t, f).*conj(J(hW(:,j), delta_t, f));
        dse_Z(j, :) = J(hZ(:,j), delta_t, f).*conj(J(hZ(:,j), delta_t, f));
        csde_WZ(j, :) = J(hW(:,j), delta_t, f).*conj(J(hZ(:,j), delta_t, f));
    end
    
    % apply multitaper
    dse1_mt = transpose(1./(1:2*N*w)).*cumsum(dse1);
    dse2_mt = transpose(1./(1:2*N*w)).*cumsum(dse2);
    csde12_mt = transpose(1./(1:2*N*w)).*cumsum(csde12);
    % apply multitaper to AR processes
    dse_X_mt = transpose(1./(1:2*N*w)).*cumsum(dse_X);
    dse_Y_mt = transpose(1./(1:2*N*w)).*cumsum(dse_Y);
    csde_XY_mt = transpose(1./(1:2*N*w)).*cumsum(csde_XY);
    % apply multitaper to dependent processes
    dse_W_mt = transpose(1./(1:2*N*w)).*cumsum(dse_W);
    dse_Z_mt = transpose(1./(1:2*N*w)).*cumsum(dse_Z);
    csde_WZ_mt = transpose(1./(1:2*N*w)).*cumsum(csde_WZ);

    % find coherence
    gamma_hat_WN = [];
    gamma_hat = [];
    gamma_hat_WZ = [];
    
    for j = (1:2*N*w)
        gamma_hat_WN(j, :) = abs(csde12_mt(j, :)).^2./(dse1_mt(j, :).*dse2_mt(j, :));
        gamma_hat(j, :) = abs(csde_XY_mt(j, :)).^2./(dse_X_mt(j, :).*dse_Y_mt(j, :));
        gamma_hat_WZ(j, :) = abs(csde_WZ_mt(j, :)).^2./(dse_W_mt(j, :).*dse_Z_mt(j, :));
    end
    
    gamma3_f_WN(i) = gamma_hat_WN(K, n);
    gamma3_f(i) = gamma_hat(K, n);
    gamma3_f_WZ(i) = gamma_hat_WZ(K, n);
end

subplot(3, 3, 1)
plot(Goodman_QQ_Plots(0, gamma1_f_WN, K), sort(gamma1_f_WN), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
ylabel('Empirical')
title('White noise')
subplot(3, 3, 2)
plot(Goodman_QQ_Plots(0, gamma1_f, K), sort(gamma1_f), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
xlabel('Theoretical, f = 3/128');
title('Independent AR(2) and AR(4)')
subplot(3, 3, 3)
plot(Goodman_QQ_Plots(gamma(n), gamma1_f_WZ, K), sort(gamma1_f_WZ), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
title('Dependent processes')

subplot(3, 3, 4)
plot(Goodman_QQ_Plots(0, gamma2_f_WN, K), sort(gamma2_f_WN), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
ylabel('Empirical')
subplot(3, 3, 5)
plot(Goodman_QQ_Plots(0, gamma2_f, K), sort(gamma2_f), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
xlabel('Theoretical, f = 50/128');
subplot(3, 3, 6)
plot(Goodman_QQ_Plots(gamma(n), gamma2_f_WZ, K), sort(gamma2_f_WZ), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off

subplot(3, 3, 7)
plot(Goodman_QQ_Plots(0, gamma3_f_WN, K), sort(gamma3_f_WN), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
ylabel('Empirical')
subplot(3, 3, 8)
plot(Goodman_QQ_Plots(0, gamma3_f, K), sort(gamma3_f), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off
xlabel('Theoretical, f = 24/1024');
subplot(3, 3, 9)
plot(Goodman_QQ_Plots(gamma(n), gamma3_f_WZ, K), sort(gamma3_f_WZ), 'p','Color','k')
hold on
plot([0 1], [0 1], 'LineWidth', 1,'Color','r');
hold off