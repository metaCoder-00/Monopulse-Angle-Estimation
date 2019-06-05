num_array = 8;
f = 10e9;
c = 3e8;
wavelength = c/f;
margin = wavelength/2;
theta_signal = [-5; 5] + 0;

dir1 = exp(-1j*2*pi*margin*(0: num_array - 1)'*sind(theta_signal(1))/wavelength);
dir2 = exp(-1j*2*pi*margin*(0: num_array - 1)'*sind(theta_signal(2))/wavelength);

theta = (-80: 0.1: 80)';
pattern1 = zeros(length(theta), 1);
pattern2 = zeros(length(theta), 1);
for n = 1: length(theta)
    steerVec = exp(-1j*2*pi*margin*(0: num_array - 1)'*sind(theta(n))/wavelength);
    pattern1(n) = dir1'*steerVec;
    pattern2(n) = dir2'*steerVec;
end
sum_beam = pattern2 + conj(pattern1);
dif_beam = pattern2 - conj(pattern1);
figure(1)
plot(theta, 20*log10(abs(sum_beam)/max(abs(sum_beam))))
hold on
plot(theta, 20*log10(abs(dif_beam)/max(abs(dif_beam))))
hold off
legend('sum beam', 'difference beam')
set(gca, 'XTICK', -75: 15: 75)
grid on
xlabel('angle/degree')
ylabel('spectrum/dB')
title('Pattern (\theta_c=0^\circ)')

ratio = dif_beam./sum_beam;
figure(2)
plot(theta(find(theta == -5): find(theta == 5)),...
    real(ratio(find(theta == -5): find(theta == 5))))
grid on
xlabel('angle/degree')
ylabel('amplitude')
title('Amplitude-Angle response')