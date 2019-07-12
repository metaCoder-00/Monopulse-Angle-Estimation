clear
SENSOR_NUM = 8;
MARGIN = 0.5;
SNR = 10;
JNR = 15;
SNAPSHOTS = 1000;
BEAM_DIR = 20;

theta_s = 25;
theta_j = 18;
amp_s = sqrt(10^(SNR/10));
% amp_j = 0;
amp_j = sqrt(10^(JNR/10));

f = 10e6;
fs = 2.5*f;
Ts = (0:SNAPSHOTS - 1)'/fs;

signal = amp_s*exp(1j*2*pi*f*Ts + 2*pi*rand(SNAPSHOTS, 1));
jammer = amp_j*exp(1j*2*pi*f*Ts + 2*pi*rand(SNAPSHOTS, 1));
noise = randn(SENSOR_NUM, SNAPSHOTS) + 1j*randn(SENSOR_NUM, SNAPSHOTS);
sv_s = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sind(theta_s));
sv_j = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sind(theta_j));
jammer_noise = sv_j*jammer.' + noise;
samples = sv_s*signal.' + jammer_noise;
covMat = jammer_noise*jammer_noise'/SNAPSHOTS;

sv_dir = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sind(BEAM_DIR));
w_sum = pinv(covMat)*sv_dir/sqrt(sv_dir'*pinv(covMat)*sv_dir);
dSv_dir = (-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)').*sv_dir;
w_dif = pinv(covMat)*dSv_dir*sqrt(sv_dir'*pinv(covMat)*sv_dir);

c_x = ((w_sum'*sv_dir)*dSv_dir - (w_sum'*dSv_dir)*sv_dir)/(w_sum'*sv_dir)^2;
B = w_dif'*c_x;
Q = covMat - (sv_dir*w_sum'*covMat)/(w_sum'*sv_dir) - (covMat*w_sum*sv_dir')/(sv_dir'*w_sum) ...
    + ((sv_dir*sv_dir')*(w_sum'*covMat*w_sum))/abs(w_sum'*sv_dir)^2;
Q1 = w_dif'*Q*w_dif;
samples = mean(samples, 2);
R = (w_dif'*samples)/(w_sum'*samples) - (w_dif'*sv_dir)/(w_sum'*sv_dir);

u_hat = real(pinv(real(B'*pinv(Q1)*B))*B'*pinv(Q1)*R);
theta_hat = BEAM_DIR + asind(u_hat)
