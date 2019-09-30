clear
array_main = 16;        % Main array
array_auxi = 4;         % Auxiliary array
array_num = array_main + array_auxi;
d = 0.5;

SNR = 15;
JNR = [20; 25];
SNAPSHOTS = 100;

boresight = 0;
theta_j = [20, -15];       % Jammer direction

sv_bore = steervec((0:array_num - 1)*d, boresight);
jammer = steervec((0:array_num - 1)*d, theta_j)*...
         sqrt(10.^(JNR/20).*exp(1j*2*pi*rand(length(theta_j), SNAPSHOTS)));
noise = randn(array_num, SNAPSHOTS) + 1j*randn(array_num, SNAPSHOTS);
data = jammer + noise;       % Data collected
data_m = data(1:array_main, :);              % Main array data
data_a = data(array_main + 1:end, :);        % Auxiliary array data
cov_mat_aa = data_a*data_a'/SNAPSHOTS;          % Auxiliary covariance matrix
cov_mat_am = data_a*data_m'/SNAPSHOTS;          % Covariance matrix of auxiliary & main array

w_s = sv_bore(1:array_main);        % Main array Sigma weight
w_d = [-ones(array_main/2, 1); ones(array_main/2, 1)].*w_s;

%-----Blocking matrix-----%
block_u = sind(-5);      % Block upper bound
block_l = sind(5);     % Block lower bound
P = zeros(array_auxi);
for m = 1:array_auxi
    for n = 1:array_auxi
        P(m, n) = (block_u - block_l)*0.5*sinc((m - n)*(block_u - block_l)/2);
    end
end
D = diag(exp(-1j*pi*(0:array_auxi - 1)'*(block_u + block_l)*0.5));
[~, S, V] = svd(P);
B = eye(array_auxi) - D*V(:, 1:array_auxi/2)*V(:, 1:array_auxi/2)'*D';      % Blocking matrix

w_s_a = pinv(B*cov_mat_aa*B')*B*cov_mat_am*w_s;

dTheta = [-4; 0; 4];        % Linear interval [-4,4]
k = cosd(boresight);        % Slope
A_a = zeros(array_auxi, length(dTheta));
A_m = zeros(array_main, length(dTheta));
rho = zeros(1, length(dTheta));
for n = 1:length(dTheta)
    sv = steervec((0:array_num - 1)*d, boresight + dTheta(n));
    A_m(:, n) = sv(1:array_main);
    A_a(:, n) = sv(array_main + 1:end);
    rho(n) = k*dTheta(n)*(w_s'*sv(1:array_main) - w_s_a'*B*sv(array_main + 1:end));
end

lag = 2*pinv(A_a'*pinv(B*cov_mat_aa*B')*A_a)*...
      (A_m'*w_d - rho' - A_a'*pinv(B*cov_mat_aa*B')*B*cov_mat_am*w_d);      % Lagrange multiplier
w_d_a = pinv(B*cov_mat_aa*B')*(B*cov_mat_am*w_d + 0);
    
% theta_s = boresight + (-4:4);
% MC_L = 1000;
% RMSE = zeros(length(theta_s), 1);
% for n = 1:length(theta_s)
%     for m = 1:MC_L
%         signal = steervec((0:array_num - 1)*d, theta_s(n))* ...
%                  sqrt(10^(SNR/20))*exp(1j*2*pi*0);
%         jammer = steervec((0:array_num - 1)*d, theta_j)* ...
%                  sqrt(10.^(JNR/20).*exp(1j*2*pi*rand(length(theta_j), SNAPSHOTS)));
%         noise = randn(array_num, SNAPSHOTS) + 1j*randn(array_num, SNAPSHOTS);
%         x = signal + jammer + noise;
%         x = mean(x, 2);
%         ratio = (w_d'*x(1:array_main) - w_d_a'*B*x(array_main + 1:end))/ ...
%                 (w_s'*x(1:array_main) - w_s_a'*B*x(array_main + 1:end));
%         theta_hat = ratio/k;
%         RMSE(n) = RMSE(n) + abs(theta_hat - theta_s(n));
%     end
%     RMSE(n) = RMSE(n)/MC_L;
% end
% figure
% plot(theta_s, RMSE)
% grid on
% xlabel('\theta (\circ)')
% ylabel('RMSE (\circ)')
% title('RMSE of JLC')

%-----DEBUG-----%
theta = (-85:0.1:85)';
pattern_m = zeros(length(theta), 1);
pattern_a = zeros(length(theta), 1);
pattern_s = zeros(length(theta), 1);
pattern_d = zeros(length(theta), 1);
for n = 1:length(theta)
    sv = steervec((0:array_num - 1)*d, theta(n));
    pattern_m(n) = w_d'*sv(1:array_main);
    pattern_a(n) = w_d_a'*B*sv(array_main + 1:end);
    pattern_s(n) = w_s'*sv(1:array_main) - w_s_a'*B*sv(array_main + 1:end);
    pattern_d(n) = w_d'*sv(1:array_main) - w_d_a'*B*sv(array_main + 1:end);
end
angle_idx = find(theta == -4):find(theta == 4);
MRC = imag(pattern_d(angle_idx)./pattern_s(angle_idx));
pattern_m = 20*log10(abs(pattern_m)/max(abs(pattern_m)));
pattern_a = 20*log10(abs(pattern_a)/max(abs(pattern_a)));
pattern_s = 20*log10(abs(pattern_s)/max(abs(pattern_s)));
pattern_d = 20*log10(abs(pattern_d)/max(abs(pattern_d)));
pattern_m = max(-50, pattern_m);
pattern_a = max(-50, pattern_a);
pattern_d = max(-50, pattern_d);
figure
plot(theta, pattern_m)
hold on 
plot(theta, pattern_a)
legend('\Delta_m', '\Delta_a')
grid on
set(gca, 'XTICK', -85: 15: 85)
axis([-90, 90, -50, 0])
xlabel('\theta (\circ)')
ylabel('Power (dB)')
title('Blocking matrix using auxiliary array')
figure
plot(theta, pattern_s)
hold on 
plot(theta, pattern_d)
grid on
legend('\Sigma', '\Delta')
set(gca, 'XTICK', -85: 5: 85)
axis([-90, 90, -50, 0])
xlabel('\theta (\circ)')
ylabel('Power (dB)')
figure
plot(theta(angle_idx), MRC)
grid on
xlabel('\theta (\circ)')
ylabel('\Delta/\Sigma')
title('Monopulse Ratio Curve')