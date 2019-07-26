num_array = 100;
num_subarray = 4;
margin = 0.5;
BEAM_DIR = 0;

NBAR = 12;
SLL = -25;

sv_beam = exp(-1j*2*pi*margin*(0:num_array - 1)'*sind(BEAM_DIR));
w_sum = taylorwin(num_array, NBAR, -35);
w_dif = DE_Dif(num_subarray, num_array/2, w_sum, SLL);

theta = (-85:0.1:85)';
sum_beam = zeros(length(theta), 1);
dif_beam = zeros(length(theta), 1);
for n = 1:length(theta)
    sv = exp(-1j*2*pi*margin*(0:num_array - 1)'*sind(theta(n)));
    sum_beam(n) = (w_sum.*sv_beam)'*sv;
    dif_beam(n) = (w_dif.*sv_beam)'*sv;
end

% plot(theta, 20*log10(abs(sum_beam)/max(abs(sum_beam))))
% grid on
% axis([-90, 90, -50, 0])
% set(gca, 'XTICK', -90: 15: 90)
% set(gca, 'YTICK', -50: 10: 0)
plot(theta, 20*log10(abs(dif_beam)/max(abs(dif_beam))))
grid on
axis([-90, 90, -50, 0])
set(gca, 'XTICK', -90: 15: 90)
set(gca, 'YTICK', -50: 10: 0)

function w = DE_Dif(P, N, w_sum, SLL)
    CR = 0.7;
    F = 0.5;
    f_th = 1e-4;
    MAX_ITR = 1000;
    POP_SIZE = 10;

    Pop = rand(P + N, POP_SIZE);
    opt_chromo = zeros(P + N, 1);
    flag = 0;
    for n = 1:MAX_ITR
        Pop = Pop(:, randperm(POP_SIZE));
        mutant = F*(Pop(:, 1) - Pop(:, 2)) + Pop(:, 3);
        for col = 1:POP_SIZE
            newChromo = zeros(P + N, 1);
            for row = 1:P
                if rand() <= CR
                    newChromo(row) = mutant(row);
                else
                    newChromo(row) = Pop(row, col);
                end
            end
            for row = P + 1:N
                if rand() <= CR
                    newChromo(row) = floor(mutant(row) + 0.5);
                else
                    newChromo(row) = floor(Pop(row, col) + 0.5);
                end
            end
            thisCost = CostF(Pop(:, col), w_sum, SLL, P);
            newCost = CostF(newChromo, w_sum, SLL, P);
            if newCost < f_th
                opt_chromo = newChromo;
                flag = 1;
                break;
            elseif newCost < thisCost
                Pop(:, col) = newChromo;
            end
        end
        if flag == 1
            break;
        end
    end
    coef_sum = w_sum(length(w_sum)/2 + 1:end);
    if flag ~= 1
        cost = zeros(POP_SIZE, 1);
        for n = 1:POP_SIZE
            cost(n) = CostF(Pop(:, n), w_sum, SLL, P);
        end
        [~, index] = min(cost);
        opt_chromo = Pop(:, index);
    end
    w = zeros(length(coef_sum), 1);
    for n = P + 1:P + N
        delta = opt_chromo(n) == (1:P)';
        w(n - P) = coef_sum(n - P)*sum(delta.*opt_chromo(1:P));
    end
    w = [-flip(w); w];
end

function loss = CostF(w_coef, w_sum, SLL, P)
    coef_sum = w_sum(length(w_sum)/2 + 1:end);
    w_dif = zeros(length(coef_sum), 1);
    for n = P + 1:length(w_coef)
        delta = w_coef(n) == (1:P)';
        w_dif(n - P) = coef_sum(n - P)*sum(delta.*w_coef(1:P));
    end
    w_dif = [-flip(w_dif); w_dif];
    theta = (0:0.1:30)';
    dif_beam = zeros(length(theta), 1);
    for n = 1:length(theta)
        sv = exp(-1j*2*pi*0.5*(0:length(w_dif) - 1)'*sind(theta(n)));
        dif_beam(n) = w_dif'*sv;
    end
    dif_beam = 20*log10(abs(dif_beam)/max(abs(dif_beam)));
    peaks = findpeaks(dif_beam);
    if isempty(peaks)
        loss = inf;
    else
        thisSLL = peaks(2);
        loss = (thisSLL - SLL)^2*heaviside(thisSLL - SLL);
    end
end