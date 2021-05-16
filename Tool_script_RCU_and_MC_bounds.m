

% This script is to plot the RCU bound and the MC bound using saddlepoint
% approximation.

set(0,'DefaultTextFontName','Times','DefaultTextFontSize',16,...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',14,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


gamma_s = 1:0.1:5;  % the definition aligns with the journal manuscript

% % CRC-ZTCC setup
% k = 64;
% m = 3;
% v = 3;
% omega = 2;
% n = omega*(k + m + v);
% R = k / n;


% CRC-TBCC setup
k = 3;
n = 4;
v = 6;
mu = ceil(v/k);
K = 144;
m = 12;
N = (K+m+k*mu)/k*n;
R = K / N;


pX = [.5; .5];
X = [-1; 1];
EbN0 = gamma_s/(2*R); % E_b/N_0 = gamma_s/(2*R)


rcu_bounds = zeros(1, size(gamma_s, 2));
mc_bounds = zeros(1, size(gamma_s, 2));


for ii = 1:size(gamma_s, 2)
    SNR = gamma_s(ii);
    sigma2 = 10^(-SNR/10);
    rcu_bounds(ii) = rcu(N, R, pX, X, sigma2);
    mc_bounds(ii) = mc(N, R, pX, X, sigma2);
end


path = './Simulation_results/';
% timestamp = datestr(now, 'mmddyy_HHMMSS');

fileName = ['RCU_and_MC_bound_N_',num2str(N),'_K_',num2str(K),'.mat'];
save([path, fileName], 'gamma_s','rcu_bounds','mc_bounds');

%%
figure;
semilogy(gamma_s, rcu_bounds, '-'); hold on
semilogy(gamma_s, mc_bounds, '-');hold on
% xline(2.5, '--k');
grid on
ylim([10^-7, 1]);
legend('RCU bound','MC bound');
xlabel('$\gamma_s$ (dB)', 'interpreter','latex');
ylabel('Probability of error');
title('K = 144, N = 216, R = 2/3');


