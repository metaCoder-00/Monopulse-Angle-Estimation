clear
M = 32;         % Main array elements number
M_auxi = 8;     % Auxiliary array

boresight = 0;
SNAPSHOTS = 100;

w_s = exp(-1j*2*pi*0.5*(0:M - 1)'*sind(boresight));

noise = randn(M + M_auxi, SNAPSHOTS) + 1j*randn(M + M_auxi, SNAPSHOTS);
x = noise;
x_m = x(1:M, :);                    % Main array data
x_a = x(M + 1:M + M_auxi, :);       % Auxiliary array data
cov_mat_a = x_a*x_a'/SNAPSHOTS;     % Covariance matrix of auxiliary array
cov_mat_am = x_a*x_m'/SNAPSHOTS;    % Covariance matrix of main array and auxiliary

%-----Form blocking matrix-----%
block_u = sind(-15);      % Block upper bound
block_l = sind(15);     % Block lower bound
P = zeros(M_auxi);
for m = 1:M_auxi
    for n = 1:M_auxi
        P(m, n) = (block_u - block_l)*0.5*sinc((m - n)*(block_u - block_l)/2);
    end
end
D = diag(exp(-1j*pi*(0:M_auxi - 1)'*(block_u + block_l)*0.5));
[~, S, V] = svd(P);
B = eye(M_auxi) - D*V(:, 1:M_auxi/2)*V(:, 1:M_auxi/2)'*D';      % Blocking matrix
w_s_a = pinv(B*cov_mat_a*B')*B*cov_mat_am*w_s;

theta = (-85:0.1:85)';
pattern_m = zeros(length(theta), 1);
pattern_a = zeros(length(theta), 1);
for n = 1:length(theta)
    sv = exp(-1j*2*pi*0.5*(0:M + M_auxi - 1)'*sind(theta(n)));
    pattern_m(n) = w_s'*sv(1:M);
    pattern_a(n) = w_s_a'*sv(M + 1:end);
end
pattern_m = 20*log10(abs(pattern_m)/max(abs(pattern_m)));
pattern_a = 20*log10(abs(pattern_a)/max(abs(pattern_a)));
pattern_m = max(-70, pattern_m);
pattern_a = max(-70, pattern_a);
figure
plot(theta, pattern_m)
hold on 
plot(theta, pattern_a)
grid on
set(gca, 'XTICK', -75: 15: 75)
axis([-85, 85, -70, 0])
xlabel('\theta (\circ)')
ylabel('Power (dB)')
title('Blocking matrix using auxiliary array')
