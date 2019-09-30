clear
array_num = 20;
subarray_num = 20;
d = 0.5;        % Margin of two adjacent elements

SNR = 15;
JNR = [20; 25];
SNAPSHOTS = 100;

boresight = 0;
theta_j = [20, -15];       % Jammer direction

T = zeros(array_num, subarray_num);
for col = 1:subarray_num
    sv = steervec((0:array_num/subarray_num - 1)*d, boresight);
    T((col - 1)*(array_num/subarray_num) + 1:col*(array_num/subarray_num), col) = sv;
end

sv_bore = T'*steervec((0:array_num - 1)*d, boresight);
jammer = steervec((0:array_num - 1)*d, theta_j)* ...
         (sqrt(10.^(JNR/20)).*exp(1j*2*pi*rand(length(theta_j), SNAPSHOTS)));
noise = randn(array_num, SNAPSHOTS) + 1j*randn(array_num, SNAPSHOTS);
data_jn = T'*(jammer + noise);       % Data without signal
cov_mat = data_jn*data_jn'/SNAPSHOTS;       % Covariance matrix of noise plus jammer
inv_cov_mat = pinv(cov_mat);

w_s = inv_cov_mat*sv_bore/(sv_bore'*inv_cov_mat*sv_bore);       % Adaptive Sigma weight

%-----Joint Linear Constraints-----%
dTheta = [-4; 0; 4];        % Linear interval [-4,4]
k = cosd(boresight);        % Slope
H = zeros(subarray_num, length(dTheta));
rho = zeros(1, length(dTheta));
for n = 1:length(dTheta)
    H(:, n) = T'*steervec((0:array_num - 1)*d, boresight + dTheta(n));
    rho(n) = k*dTheta(n)*w_s'*T'*steervec((0:array_num - 1)*d, boresight + dTheta(n));
end
w_d = inv_cov_mat*H*pinv(H'*inv_cov_mat*H)*rho';        % Adaptive Delta weight

theta_s = boresight + (-4:4)';      % Signal direction
MC_L = 1000;      
RMSE = zeros(length(theta_s), 1);
for n = 1:length(theta_s)
    for m = 1:MC_L
        signal = steervec((0:array_num - 1)*d, theta_s(n))* ...
                 sqrt(10^(SNR/20))*exp(1j*2*pi*rand(length(theta_s(n)), SNAPSHOTS));
        jammer = steervec((0:array_num - 1)*d, theta_j)* ...
                 sqrt(10.^(JNR/20).*exp(1j*2*pi*rand(length(theta_j), SNAPSHOTS)));
        noise = randn(array_num, SNAPSHOTS) + 1j*randn(array_num, SNAPSHOTS);
        x = T'*(signal + jammer + noise);
        x = mean(x, 2);
        ratio = real((w_d'*x)/(w_s'*x));
        theta_hat = ratio/k;
        RMSE(n) = RMSE(n) + abs(theta_hat - theta_s(n));
    end
    RMSE(n) = RMSE(n)/MC_L;
end
figure
plot(theta_s, RMSE)
grid on
xlabel('\theta (\circ)')
ylabel('RMSE (\circ)')
title('RMSE of JLC')

%%
%-----DEBUG-----%
theta = (-85:0.1:85)';
pattern_s = zeros(length(theta), 1);
pattern_d = zeros(length(theta), 1);
for n = 1:length(theta)
    sv = T'*steervec((0:array_num - 1)*d, theta(n));
    pattern_s(n) = w_s'*sv;
    pattern_d(n) = w_d'*sv;
end
angle_idx = find(theta == -4):find(theta == 4);
MRC = real(pattern_d(angle_idx)./pattern_s(angle_idx));
pattern_s = 20*log10(abs(pattern_s)/max(abs(pattern_s)));
pattern_d = 20*log10(abs(pattern_d)/max(abs(pattern_d)));
pattern_s = max(-50, pattern_s);
pattern_d = max(-50, pattern_d);
figure
plot(theta, pattern_s)
hold on
plot(theta, pattern_d)
grid on
legend('\Sigma', '\Delta')
xlabel('\theta (\circ)')
ylabel('Power (dB)')
title('Adaptive \Sigma & \Delta Pattern')
figure
plot(theta(angle_idx), MRC)
grid on
xlabel('\theta (\circ)')
ylabel('\Delta/\Sigma')
title('Monopulse Ratio Curve')