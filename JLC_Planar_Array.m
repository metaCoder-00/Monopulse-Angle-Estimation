ROW_NUM = 16;
COL_NUM = 16;
ROW_PART = 2;
COL_PART = 2;
SUB_ROW_NUM = ROW_NUM/ROW_PART;
SUB_COL_NUM = COL_NUM/COL_PART;
SIR = 5;
SNR = 10;
SNAPSHOTS = 100;

f = 10e9;
c = 3e8;
wavelength = c/f;
margin_row = wavelength/2;
margin_col = wavelength/2;
direction = [20; 20];
jammer_main = direction + [2; -2];
jammer_side = direction + [-50; -50];

transMat_u = zeros(COL_NUM, COL_PART);
for col = 1:COL_PART
    steerVec = exp(-1j*2*pi*margin_col*(0:SUB_COL_NUM - 1)'* ... 
               sind(direction(1))*cosd(direction(2))/wavelength);
    transMat_u((col - 1)*SUB_COL_NUM + 1:col*SUB_COL_NUM, col) = steerVec;
end
transMat_v = zeros(ROW_NUM, ROW_PART);
for col = 1:ROW_PART
    steerVec = exp(-1j*2*pi*margin_row*(0:SUB_ROW_NUM - 1)'* ... 
               sind(direction(2))/wavelength);
    transMat_v((col - 1)*SUB_ROW_NUM + 1:col*SUB_ROW_NUM, col) = steerVec;
end
transMat = kron(transMat_u, transMat_v);

steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, direction(1), direction(2));
samples = zeros(ROW_NUM*COL_NUM, SNAPSHOTS);
for n = 1:SNAPSHOTS
    jammerMainVec = zeros(ROW_NUM, COL_NUM);
    jammerSideVec = zeros(ROW_NUM, COL_NUM);
    noise = zeros(ROW_NUM, COL_NUM);
    for row = 1:ROW_NUM
        for col = 1:COL_NUM        
            jammer_main_phase_shift = 2*pi*((row - 1)*margin_row* ... 
                  sind(jammer_main(1))*cosd(jammer_main(2)) + ...
                  (col - 1)*margin_col*sind(jammer_main(2)))/wavelength;
            jammerMainVec(row, col) = sqrt(1/10^(SIR/10))*randn()*exp(-1j*jammer_main_phase_shift);
            jammer_side_phase_shift = 2*pi*((row - 1)*margin_row* ... 
                                      sind(jammer_side(1))*cosd(jammer_side(2)) + ...
                                      (col - 1)*margin_col*sind(jammer_side(2)))/wavelength;
            jammerSideVec(row, col) = sqrt(1/10^(SIR/10))*randn()*exp(-1j*jammer_side_phase_shift);
            noise(row, col) = sqrt(1/10^(SNR/10))*randn();
        end
    end
    jammerMainVec = jammerMainVec(:);
    jammerSideVec = jammerSideVec(:);
    noise = noise(:);
    data = steerVec + jammerMainVec + jammerSideVec + noise;
    samples(:, n) = data;
end

sub_samples = transMat'*samples;
weight_sum = transMat'*steerVec;
coef = [
        pi*ROW_NUM*margin_row*cosd(direction(2))
        pi*COL_NUM*margin_col*cosd(direction(1))*cosd(direction(2))
        -pi*COL_NUM*margin_col*sind(direction(1))*sind(direction(2))
        ]/wavelength;
covMat = sub_samples*sub_samples'/SNAPSHOTS;
weight_sum = pinv(covMat)*weight_sum;
delta = [3; 3];

sub_sv = transMat'*steerVec;

azm = direction(1) + delta(1);
elv = direction(2) + delta(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_ap_ep = transMat'*steerVec;

azm = direction(1) + delta(1);
elv = direction(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_ap = transMat'*steerVec;

azm = direction(1) + delta(1);
elv = direction(2) - delta(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_ap_em = transMat'*steerVec;

azm = direction(1);
elv = direction(2) + delta(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_ep = transMat'*steerVec;

azm = direction(1);
elv = direction(2) - delta(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_em = transMat'*steerVec;

azm = direction(1) - delta(1);
elv = direction(2) + delta(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_am_ep = transMat'*steerVec;

azm = direction(1) - delta(1);
elv = direction(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_am = transMat'*steerVec;

azm = direction(1) - delta(1);
elv = direction(2) - delta(2);
steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv);
sub_sv_am_em = transMat'*steerVec;

consMat = [sub_sv_ap_ep, sub_sv_ap, sub_sv_ap_em, ...
           sub_sv_ep, sub_sv, sub_sv_em, ...
           sub_sv_am_ep, sub_sv_am, sub_sv_am_em];
consVec_a = [(coef(2)*delta(1) + coef(3)*delta(2))*weight_sum'*sub_sv_ap_ep, ...
             (coef(2)*delta(1))*weight_sum'*sub_sv_ap, ...
             (coef(2)*delta(1) - coef(3)*delta(2))*weight_sum'*sub_sv_ap_em, ...
             (coef(3)*delta(2))*weight_sum'*sub_sv_ep, ...
             0, ...
             (-coef(3)*delta(2))*weight_sum'*sub_sv_em, ...
             (-coef(2)*delta(1) + coef(3)*delta(2))*weight_sum'*sub_sv_am_em, ...
             (-coef(2)*delta(1))*weight_sum'*sub_sv_am, ...
             (-coef(2)*delta(1) - coef(3)*delta(2))*weight_sum'*sub_sv_am_em];
consVec_e = [(coef(1)*delta(2))*weight_sum'*sub_sv_ap_ep, ...
             0, ...
             -(coef(1)*delta(2))*weight_sum'*sub_sv_ap_em, ...
             (coef(1)*delta(2))*weight_sum'*sub_sv_ep, ...
             0, ...
             -(coef(1)*delta(2))*weight_sum'*sub_sv_em, ...
             (coef(1)*delta(2))*weight_sum'*sub_sv_am_ep, ...
             0, ...
             -(coef(1)*delta(2))*weight_sum'*sub_sv_am_em];

weight_dif_a = pinv(covMat)*consMat*pinv(consMat'*pinv(covMat)*consMat)*consVec_a';
weight_dif_e =  pinv(covMat)*consMat*pinv(consMat'*pinv(covMat)*consMat)*consVec_e';

theta = (17:0.1:23)';
phi = (17:0.1:23)';
sum_beam = zeros(length(theta), length(phi));
dif_beam_a = zeros(length(theta), length(phi));
dif_beam_e = zeros(length(theta), length(phi));
for m = 1:length(theta)
    for n = 1:length(phi)
        steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, theta(m), phi(n));
        sub_sv = transMat'*steerVec;
        sum_beam(m, n ) = weight_sum'*sub_sv;
        dif_beam_a(m, n) = weight_dif_a'*sub_sv;
        dif_beam_e(m, n) = weight_dif_e'*sub_sv;
    end
end

figure
plot3(theta, phi, imag(dif_beam_a./sum_beam))
figure
plot3(theta, phi, imag(dif_beam_e./sum_beam))

function steerVec = planarSteerVec(ROW_NUM, COL_NUM, margin_row, margin_col, wavelength, azm, elv)
    steerVec = zeros(ROW_NUM, COL_NUM);
    for row = 1:ROW_NUM
        for col = 1:COL_NUM
            phase_shift = 2*pi*((row - 1)*margin_row* ... 
                          sind(azm)*cosd(elv) + ...
                          (col - 1)*margin_col*sind(elv))/wavelength;
            steerVec(row, col) = exp(-1j*phase_shift);
        end
    end
    steerVec = steerVec(:);
end
