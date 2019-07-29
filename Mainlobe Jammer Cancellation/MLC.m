%
%   References:
%     [1] Kai-Bor Yu and D. J. Murrow. 
%         Adaptive digital beamforming for angle estimation in jamming.
%         IEEE Transactions on Aerospace and Electronic Systems
%         vol. 37, no. 2, pp. 508-523, April 2001.

clear
M = 16;     % Array elements number
N = 16;     % Array elements number
SNR = 15;
MLJNR = 20; % Mainlobe jammer to noise ratio
SLJNR = 20; % Sidelobe jammer to noise ratio
SNAPSHOTS = 100;

boresight = [0; 45];
mainlobe_jammer = [3; 42];
sidelobe_jammer = [30, 15];

sv_dir = planar_steervec(M, N, boresight);                  % Boresight steer vector 
w_Sa = kron(taylorwin(M, M/2, -35), ones(N, 1)).*sv_dir;   % Taylor weight @ azimuth
w_Se = kron(ones(M, 1), taylorwin(N, N/2, -35)).*sv_dir;   % Taylor weight @ elevation
w_Da = kron(baylisswin(M, M/2, -35), ones(N, 1)).*sv_dir;  % Bayliss weight @ azimuth
w_De = kron(ones(M, 1), baylisswin(N, N/2, -25)).*sv_dir;  % Bayliss weight @ elevation

%----------Generate samples--------%
MLJ = sqrt(10^(MLJNR/20))*planar_steervec(M, N, mainlobe_jammer);
SLJ = sqrt(10^(SLJNR/20))*planar_steervec(M, N, sidelobe_jammer);
noise = randn(M*N, SNAPSHOTS) + 1j*randn(M*N, SNAPSHOTS);
samples = MLJ + SLJ + noise;

data_S = (w_Sa'*samples).*(w_Se'*samples);        % Sum pattern
data_DA = (w_Da'*samples).*(w_Se'*samples);      % Difference pattern of azimuth
data_DE = (w_Sa'*samples).*(w_De'*samples);      % Difference pattern of elevation
data_DD = (w_Da'*samples).*(w_De'*samples);      % Double difference pattern

corrMat_S_DA = data_S*data_DA'/SNAPSHOTS;       % Correlation matrix of sum and Delta-azimuth pattern
corrMat_DA_DA = data_DA*data_DA'/SNAPSHOTS;     % Correlation matrix of Delta-azimuth and Delta-azimuth pattern
corrMat_DE_DD = data_DE*data_DD'/SNAPSHOTS;     % Correlation matrix of Delta-elevation and Delta-Delta pattern
corrMat_DD_DD = data_DD*data_DD'/SNAPSHOTS;     % Correlation matrix of Delta-Delta and Delta-Delta pattern

w_a = (corrMat_S_DA*pinv(corrMat_DA_DA) + corrMat_DE_DD*pinv(corrMat_DD_DD))/2;     % Mainlobe jammer canceller weight @ azimuth

corrMat_S_DE = data_S*data_DE'/SNAPSHOTS;       % Correlation matrix of sum and Delta-elevation pattern
corrMat_DE_DE = data_DE*data_DE'/SNAPSHOTS;     % Correlation matrix of Delta-elevation and Delta-elevation pattern
corrMat_DA_DD = data_DA*data_DD'/SNAPSHOTS;     % Correlation matrix of Delta-azimuth and Delta-Delta pattern

w_e = (corrMat_S_DE*pinv(corrMat_DE_DE) + corrMat_DA_DD*pinv(corrMat_DD_DD))/2;     % Mainlobe jammer canceller weight @ elevation

%%
%---------Beam ratio----------%
theta = boresight(1) + (-3:0.1:3)';
phi = boresight(2) + (-3:0.1:3)';
beam_S_a = zeros(length(theta), length(phi));
beam_S_e = zeros(length(theta), length(phi));
beam_D_a = zeros(length(theta), length(phi));
beam_D_e = zeros(length(theta), length(phi));
for m = 1:length(theta)
    for n = 1:length(phi)
        x = sqrt(10^(SNR/20))*planar_steervec(M, N, [theta(m); phi(n)]) + ...
            sqrt(10^(MLJNR/20))*planar_steervec(M, N, mainlobe_jammer) + ...
            sqrt(10^(SLJNR/20))*planar_steervec(M, N, sidelobe_jammer);
        beam_S_a(m, n) = w_Sa'*x;
        beam_S_e(m, n) = w_Se'*x;
        beam_D_a(m, n) = w_Da'*x;
        beam_D_e(m, n) = w_De'*x;
    end
end
ratio_A = (beam_D_a.*beam_S_e - w_e.*beam_D_a.*beam_D_e)./ ...
          (beam_S_a.*beam_S_e - w_e.*beam_S_a.*beam_D_e);
ratio_E = (beam_S_a.*beam_D_e - w_a.*beam_D_a.*beam_D_e)./ ...
          (beam_S_a.*beam_S_e - w_a.*beam_D_a.*beam_S_e);

figure
[x, y] = meshgrid(theta, phi);
mesh(x, y, imag(ratio_A).')
title('\Delta_A to \Sigma_A ratio')
xlabel('Azimuth/degree')
ylabel('Elevation/degree')
zlabel('\Delta_A/\Sigma_A')
figure
[x, y] = meshgrid(theta, phi);
mesh(x, y, imag(ratio_E).')
title('\Delta_E to \Sigma_E ratio')
xlabel('Azimuth/degree')
ylabel('Elevation/degree')
zlabel('\Delta_E/\Sigma_E')

function sv = planar_steervec(M, N, dir)
    sv_amz = exp(-1j*2*pi*0.5*(0:M - 1)'*sind(dir(1))*cosd(dir(2)));
    sv_elv = exp(-1j*2*pi*0.5*(0:N - 1)'*sind(dir(2)));
    sv = kron(sv_amz, sv_elv);
end