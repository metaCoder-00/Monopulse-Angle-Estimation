clear
SENSOR_NUM = 8;
MARGIN = 0.5;
SNR = 10;
JNR = 15;
SNAPSHOTS = 500;
BEAM_DIR = 20;

theta_s = 25;
theta_j = 15;
amp_s = sqrt(10^(SNR/10));
amp_j = sqrt(10^(JNR/10));
covMat_n = eye(SENSOR_NUM);

f = 10e6;
fs = 2.5*f;
Ts = (0:SNAPSHOTS - 1)'/fs;

signal = amp_s*exp(1j*2*pi*f*Ts + 2*pi*rand(SNAPSHOTS, 1));
jammer = amp_j*exp(1j*2*pi*f*Ts + 2*pi*rand(SNAPSHOTS, 1));
noise = randn(SENSOR_NUM, SNAPSHOTS) + 1j*randn(SENSOR_NUM, SNAPSHOTS);
sv_s = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sind(theta_s));
sv_j = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sind(theta_j));
samples = sv_s*signal.' + sv_j*jammer.' + noise;
covMat_n = covMat_n + (sv_j*jammer.')*(sv_j*jammer.')'/SNAPSHOTS;

ITR = size(samples, 2);
theta = BEAM_DIR;
sine = sind(theta);
theta_hat = theta;
for n = 1:ITR
    sv = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sine);
    w = pinv(sqrt(sv'*pinv(covMat_n)*sv))*pinv(covMat_n)*sv;
    dSv = (-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)').*sv;
    d2Sv = pinv(covMat_n)*dSv*sqrt(sv'*pinv(covMat_n)*sv);
    mu = real((dSv'*pinv(covMat_n)*sv)/(sv'*pinv(covMat_n)*sv));
    dLf = 2*(real((d2Sv'*samples(:, n))/(w'*samples(:, n))) - mu);
    d2Lf = 2*mu^2 - (2*d2Sv'*dSv)/(w'*sv);
    sine = sine - pinv(d2Lf)*dLf;
    theta_hat = [theta_hat; asind(abs(sine))];
end

plot((0:ITR)', theta_hat)
grid on
xlabel('Itrations')
ylabel('\theta (\circ)')
title('Singal + Jammer + Noise')
