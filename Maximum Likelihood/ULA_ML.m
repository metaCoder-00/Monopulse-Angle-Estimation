%
%   References:
%     [1] U. Nickel 
%         Monopulse estimation with adaptive arrays 
%         IEE Proceedings F - Radar and Signal Processing
%         vol. 140, no. 5, pp. 303-308, Oct. 1993.
clear
SENSOR_NUM = 8;
MARGIN = 0.5;
SNR = 15;
JNR = 20;
SNAPSHOTS = 600;
BEAM_DIR = 20;

theta_s = 22;
theta_j = 25;
amp_s = sqrt(10^(SNR/10));
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
covMat_n = jammer_noise*jammer_noise'/SNAPSHOTS;

%------MBGD-------%
BATCH = 200;
BATCH_SIZE = SNAPSHOTS/BATCH;
theta = BEAM_DIR;
sine = sind(theta);
theta_hat = theta;
for batch = 1:BATCH
    dir = 0;
    for n = 1:BATCH_SIZE
        sv = exp(-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)'*sine);
        w = pinv(sqrt(sv'*pinv(covMat_n)*sv))*pinv(covMat_n)*sv;
        dSv = (-1j*2*pi*MARGIN*(0:SENSOR_NUM - 1)').*sv;
        d2Sv = pinv(covMat_n)*dSv/sqrt(sv'*pinv(covMat_n)*sv);
        mu = real((dSv'*pinv(covMat_n)*sv)/(sv'*pinv(covMat_n)*sv));
        dLf = 2*(real((d2Sv'*samples(:, (batch - 1)*BATCH_SIZE + n))/(w'*samples(:, (batch - 1)*BATCH_SIZE + n))) - mu);
        d2Lf = 2*mu^2 - (2*d2Sv'*dSv)/(w'*sv);
        dir = dir + pinv(d2Lf)*dLf;
    end
    sine = sine - dir/BATCH_SIZE;
    theta_hat = [theta_hat; asind(abs(sine))];
end

plot((0:BATCH)', theta_hat)
grid on
xlabel('Batch')
ylabel('\theta (\circ)')
title('Singal + Jammer + Noise (jammer = 25\circ)')
