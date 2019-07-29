%
%   References:
%     [1] Z. Cheng, Z. He, X. Duan, X. Zhang and B. Liao
%         Adaptive Monopulse Approach With Joint Linear Constraints for Planar Array at Subarray Level
%         IEEE Transactions on Aerospace and Electronic Systems
%         vol. 54, no. 3, pp. 1432-1441, June 2018.
clear
M = 16;
N = 16;
L1 = 4;
L2 = 4;

dy = 0.5;
dz = 0.5;
SIR = -20;
SNR = 20;
SNAPSHOTS = 100;
amz = 20;
elv = 20;
amz_j = 22;
elv_j = 18;

jammer = planarSteer(M, N, amz_j, elv_j)*(sqrt(1/10^(SIR/10))*ones(1, SNAPSHOTS));
noise = sqrt(1/10^(SNR/10))*(randn(M*N, SNAPSHOTS) + 1j*randn(M*N, SNAPSHOTS));
samples = jammer + noise;


T_u = zeros(M, L1);
for col = 1:L1
    sv = exp(-1j*2*pi*dy*(0:M/L1 - 1)'*sind(amz)*cosd(elv));
    T_u((col - 1)*(M/L1) + 1:col*(M/L1), col) = sv;
end
T_v = zeros(N, L2);
for col = 1:L2
    sv = exp(-1j*2*pi*dz*(0:N/L2 - 1)'*sind(elv));
    T_v((col - 1)*(N/L2) + 1:col*(N/L2), col) = sv;
end
T = kron(T_u, T_v);

sub_samples = T'*samples;
covMat = sub_samples*sub_samples'/SNAPSHOTS;
sv = T'*planarSteer(M, N, amz, elv);
w_sum = (pinv(covMat)*sv)/(sv'*pinv(covMat)*sv);

% k = [pi*M*dy*cosd(elv); pi*N*dz*cosd(amz)*cosd(elv); -pi*N*dz*sind(amz)*sind(elv)];
k = [cosd(elv); cosd(amz)*cosd(elv); -sind(amz)*sind(elv)];

dAmz = [-3; 0; 3];
dElv = [-3; 0; 3];
H = zeros(L1*L2, length(dAmz)*length(dElv));
rho_a = zeros(1, length(dAmz)*length(dElv));
rho_e = zeros(1, length(dAmz)*length(dElv));
for m = 1:length(dAmz)
    for n = 1:length(dElv)
        H(:, (m - 1)*length(dElv) + n) = T'*planarSteer(M, N, amz + dAmz(m), elv + dElv(n));
        rho_a((m - 1)*length(dElv) + n) = (k(2)*dAmz(m) + k(3)*dElv(n))*w_sum'*H(:, (m - 1)*3 + n);
        rho_e((m - 1)*length(dElv) + n) = k(1)*dElv(n)*w_sum'*H(:, (m - 1)*length(dElv) + n);
    end
end

w_dif_a = pinv(covMat)*H*pinv(H'*pinv(covMat)*H)*rho_a';
w_dif_e = pinv(covMat)*H*pinv(H'*pinv(covMat)*H)*rho_e';

theta = (min(dAmz):0.1:max(dAmz))' + amz;
phi = (min(dElv):0.1:max(dElv))' + elv;
sum_beam = zeros(length(theta), length(phi));
dif_beam_a = zeros(length(theta), length(phi));
dif_beam_e = zeros(length(theta), length(phi));
for m = 1:length(theta)
    for n = 1:length(phi)
        sv = T'*planarSteer(M, N, theta(m), phi(n));
        sum_beam(m, n) = w_sum'*sv;
        dif_beam_a(m, n) = w_dif_a'*sv;
        dif_beam_e(m, n) = w_dif_e'*sv;
    end
end

[x, y] = meshgrid(theta, phi);
figure
mesh(x, y, real(dif_beam_a./sum_beam).')
% axis([min(theta), max(theta), min(phi), max(phi), -100, 100])
xlabel('Azimuth/degree')
ylabel('Elevation/degree')
zlabel('\Delta_a/\Sigma')
title('Azimuth dif-sum ratio')
figure
mesh(x, y, real(dif_beam_e./sum_beam).')
% axis([min(theta), max(theta), min(phi), max(phi), -100, 100])
xlabel('Azimuth/degree')
ylabel('Elevation/degree')
zlabel('\Delta_e/\Sigma')
title('Elevation dif-sum ratio')

function sv = planarSteer(M, N, amz, elv)
    sv_a = exp(-1j*2*pi*0.5*(0:M - 1)'*sind(amz)*cosd(elv));
    sv_e = exp(-1j*2*pi*0.5*(0:N - 1)'*sind(elv));
    sv = sv_e*sv_a.';
    sv = sv(:);
end
