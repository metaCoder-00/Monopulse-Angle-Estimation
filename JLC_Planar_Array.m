ROW_NUM = 16;
COL_NUM = 16;
ROW_PART = 2;
COL_PART = 2;
SUB_ROW_NUM = ROW_NUM/ROW_PART;
SUB_COL_NUM = COL_NUM/COL_PART;
SIR = 5;
SNR = 10;

f = 10e9;
c = 3e8;
wavelength = c/f;
margin_row = wavelength/2;
margin_col = wavelength/2;
theta_signal = [20; 20];
jammer_main = theta_signal + [2; -2];
jammer_side = theta_signal + [-50; -50];

transMat_u = zeros(COL_NUM, COL_PART);
for col = 1:COL_PART
    steerVec = exp(-1j*2*pi*margin_col*(0:SUB_COL_NUM - 1)'* ... 
               sind(theta_signal(1))*cosd(theta_signal(2))/wavelength);
    transMat_u((col - 1)*SUB_COL_NUM + 1:col*SUB_COL_NUM, col) = steerVec;
end
transMat_v = zeros(ROW_NUM, ROW_PART);
for col = 1:ROW_PART
    steerVec = exp(-1j*2*pi*margin_row*(0:SUB_ROW_NUM - 1)'* ... 
               sind(theta_signal(2))/wavelength);
    transMat_u((col - 1)*SUB_ROW_NUM + 1:col*SUB_ROW_NUM, col) = steerVec;
end
transMat = kron(transMat_u, transMat_v);

steerVec = zeros(ROW_NUM, COL_NUM);
jammerMainVec = zeros(ROW_NUM, COL_NUM);
jammerSideVec = zeros(ROW_NUM, COL_NUM);
noise = zeros(ROW_NUM, COL_NUM);
for row = 1:ROW_NUM
    for col = 1:COL_NUM
        phase_shift = 2*pi*((row - 1)*margin_row* ... 
                      sind(theta_signal(1))*cosd(theta_signal(2)) + ...
                      (col - 1)*margin_col*sind(theta_signal(2)))/wavelength;
        steerVec(row, col) = exp(-1j*phase_shift);
        
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
steerVec = steerVec(:);
steerVec_signal = steerVec;
jammerMainVec = jammerMainVec(:);
jammerSideVec = jammerSideVec(:);
noise = noise(:);
data = steerVec + jammerMainVec + jammerSideVec + noise;
