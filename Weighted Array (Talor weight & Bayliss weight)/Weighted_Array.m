clear
num_array = 8;
f = 10e9;
c = 3e8;
wavelength = c/f;
margin = wavelength/2;
theta_signal = 0;

dir = exp(-1j*2*pi*margin*(0: num_array - 1)'*sind(theta_signal)/wavelength);
w_s = taylorwin(num_array, 4, -35).*dir;
w_d = baylisswin(num_array, 4, -35).*dir;

theta = (-80: 0.1: 80)';
sum_beam = zeros(length(theta), 1);
dif_beam = zeros(length(theta), 1);
for n = 1: length(theta)
    steerVec = exp(-1j*2*pi*margin*(0: num_array - 1)'*sind(theta(n))/wavelength);
    sum_beam(n) = w_s'*steerVec;
    dif_beam(n) = w_d'*steerVec;
end
figure
plot(theta, 20*log10(abs(sum_beam)/max(abs(sum_beam))))
hold on
plot(theta, 20*log10(abs(dif_beam)/max(abs(dif_beam))))
hold off
legend('sum beam', 'difference beam')
axis([-80, 80, -50, 0])
set(gca, 'XTICK', -75: 15: 75)
grid on
xlabel('angle/degree')
ylabel('spectrum/dB')
title('Pattern (\theta_c=0^\circ)')

ratio = dif_beam./sum_beam;
figure(2)
plot(theta(find(theta == -5): find(theta == 5)),...
    imag(ratio(find(theta == -5): find(theta == 5))))
grid on
xlabel('angle/degree')
ylabel('amplitude')
title('Amplitude-Angle response')