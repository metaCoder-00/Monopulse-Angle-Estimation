num_array = 2*4;
f = 10e9;
c = 3e8;
wavelength = c/f;
margin = wavelength/2;
theta_signal = 0;

theta = (-90: 0.1: 90)';
sum_beam = zeros(length(theta), 1);
for n = 1: length(theta)
    steerVec = exp(-1j*2*pi*margin*(0: num_array - 1)'*sind(theta(n))/wavelength);
    weight = taylorwin(num_array);
    sum_beam(n) = weight'*steerVec;
end

plot(theta, 20*log10(abs(sum_beam)/max(abs(sum_beam))))
set(gca, 'XTICK', -80: 20: 80)
grid on
xlabel('angle/degree')
ylabel('spectrum/dB')
title('Pattern (\theta_c=0^\circ)')