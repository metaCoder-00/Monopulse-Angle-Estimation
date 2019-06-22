M = 16;
N = 16;
L1 = 4;
L2 = 4;

f = 10e9;
c = 3e8;
dy = 0.5;
dz = 0.5;
SIR = -20;
SNR = 20;
SNAPSHOTS = 100;
amz = 20;
elv = 20;
amz_j = 22;
elv_j = 18;

jammer = sqrt(1/10^(SIR/10))*planarSteer(M, N, amz_j, elv_j);
samples = jammer + sqrt(1/10^(SNR/10))*randn(M*N, SNAPSHOTS);

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
w_sum = pinv(covMat)*sv;

%--------Debug---------%
% w_sum = sv;
% w = kron([ones(L1/2, 1); -ones(L1/2, 1)], ones(L2, 1));
% w_dif_a = w.*w_sum;
% w = kron(ones(L1, 1), [ones(L2/2, 1); -ones(L2/2, 1)]);
% w_dif_e = w.*w_sum;
% theta = (amz - 3:0.1:amz + 3)';
% phi = (elv - 3:0.1:elv + 3);
% sum_beam = zeros(length(theta), length(phi));
% dif_beam_a = zeros(length(theta), length(phi));
% dif_beam_e = zeros(length(theta), length(phi));
% for m = 1:length(theta)
%     for n = 1:length(phi)
%         steerVec = planarSteer(M, N, theta(m), phi(n));
%         sub_sv = T'*steerVec;
%         sum_beam(m, n) = w_sum'*sub_sv;
%         dif_beam_a(m, n) = w_dif_a'*sub_sv;
%         dif_beam_e(m, n) = w_dif_e'*sub_sv;
%     end
% end
% 
% [x, y] = meshgrid(theta, phi);
% figure
% mesh(x, y, 20*log10(abs(sum_beam)/max(max(abs(sum_beam)))))
% mesh(x, y, 20*log10(abs(dif_beam_a)/max(max(abs(dif_beam_a)))))
% mesh(x, y, 20*log10(abs(dif_beam_e)/max(max(abs(dif_beam_e)))))
% mesh(x, y, imag(dif_beam_a./sum_beam))
% view(90, 0)
% figure
% mesh(x, y, imag(dif_beam_e./sum_beam))
% view(0, 0)

k = [pi*M*dy*cosd(elv); pi*N*dz*cosd(amz)*cosd(elv); -pi*N*dz*sind(amz)*sind(elv)];
dAmz = 3;
dElv = 3;
sv_ap = T'*planarSteer(M, N, amz + dAmz, elv);
sv_ep = T'*planarSteer(M, N, amz, elv + dElv);
sv_am = T'*planarSteer(M, N, amz - dAmz, elv);
sv_em = T'*planarSteer(M, N, amz, elv - dElv);
sv_ap_ep = T'*planarSteer(M, N, amz + dAmz, elv + dElv);
sv_ap_em = T'*planarSteer(M, N, amz + dAmz, elv - dElv);
sv_am_ep = T'*planarSteer(M, N, amz - dAmz, elv + dElv);
sv_am_em = T'*planarSteer(M, N, amz - dAmz, elv - dElv);
H_p = [sv_ap_ep, sv_ap, sv_ap_em];
H_0 = [sv_ep, sv, sv_em];
H_m = [sv_am_ep, sv_am, sv_am_em];
H = [H_p, H_0, H_m];
rho_a_p = [(k(2)*dAmz + k(3)*dElv)*w_sum'*sv_ap_ep, ...
         (k(2)*dAmz)*w_sum'*sv_ap, ...
         (k(2)*dAmz - k(3)*dElv)*w_sum'*sv_ap_em];
rho_a_0 = [(k(3)*dElv)*w_sum'*sv_ep, 0, (-k(3)*dElv)*w_sum'*sv_em];
rho_a_m = [(-k(2)*dAmz + k(3)*dElv)*w_sum'*sv_am_ep, ...
         (-k(2)*dAmz)*w_sum'*sv_am, ...
         (-k(2)*dAmz - k(3)*dElv)*w_sum'*sv_am_em];
rho_a = [rho_a_p, rho_a_0, rho_a_m];
rho_e_p = [(k(1)*dElv)*w_sum'*sv_ap_ep, 0, -(k(1)*dElv)*w_sum'*sv_ap_em];
rho_e_0 = [(k(1)*dElv)*w_sum'*sv_ep, 0, -(k(1)*dElv)*w_sum'*sv_em];
rho_e_m = [(k(1)*dElv)*w_sum'*sv_am_ep, 0, -(k(1)*dElv)*w_sum'*sv_am_em];
rho_e = [rho_e_p, rho_e_0, rho_e_m];

w_dif_a = pinv(covMat)*H*pinv(H'*pinv(covMat)*H)*rho_a';
w_dif_e = pinv(covMat)*H*pinv(H'*pinv(covMat)*H)*rho_e';

w_sum = w_sum/max(abs(w_sum));
w_dif_a = w_dif_a/max(abs(w_dif_a));
w_dif_e = w_dif_e/max(abs(w_dif_e));

theta = (amz - dAmz:0.1:amz + dAmz)';
phi = (elv - dElv:0.1:elv + dElv);
sum_beam = zeros(length(theta), length(phi));
dif_beam_a = zeros(length(theta), length(phi));
dif_beam_e = zeros(length(theta), length(phi));
for m = 1:length(theta)
    for n = 1:length(phi)
        s = planarSteer(M, N, theta(m), phi(n));
        sv = T'*(s + jammer + sqrt(1/10^(SNR/10))*randn(M*N, 1));
        sum_beam(m, n) = w_sum'*sv;
        dif_beam_a(m, n) = w_dif_a'*sv;
        dif_beam_e(m, n) = w_dif_e'*sv;
    end
end

[x, y] = meshgrid(theta, phi);
figure

%----------Debug-----------%
% mesh(x, y, 20*log10(abs(sum_beam)/max(max(abs(sum_beam)))))
% mesh(x, y, 20*log10(abs(dif_beam_a)/max(max(abs(dif_beam_a)))))
% mesh(x, y, 20*log10(abs(dif_beam_e)/max(max(abs(dif_beam_e)))))

mesh(x, y, imag(dif_beam_a./sum_beam))
% axis([amz - dAmz, amz + dAmz, elv - dElv, elv + dElv, -10, 10])
xlabel('Elevation/degree')
ylabel('Azimuth/degree')
view(90, 0)
figure
mesh(x, y, imag(dif_beam_e./sum_beam))
% axis([amz - dAmz, amz + dAmz, elv - dElv, elv + dElv, -10, 10])
xlabel('Elevation/degree')
ylabel('Azimuth/degree')
view(0, 0)


function sv = planarSteer(M, N, amz, elv)
    sv = zeros(N, M);
    for m = 1:M
        for n = 1:N
            phi = 2*pi*0.5*((m - 1)*sind(amz)*cosd(elv) + (n - 1)*sind(elv));
            sv(n, m) = exp(-1j*phi);
        end
    end
    sv = sv(:);
end
