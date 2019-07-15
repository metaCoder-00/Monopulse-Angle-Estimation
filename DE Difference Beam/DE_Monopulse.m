num_array = 100;
num_subarray = 4;
margin = 0.5;
BEAM_DIR = 0;

NBAR = 12;
SLL = -35;

sv_beam = exp(-1j*2*pi*margin*(0:num_array - 1)'*sind(BEAM_DIR));
w_sum = taylorwin(num_array, NBAR, SLL);
w_dif_coef = DE_Dif(num_subarray, num_array/2, w_sum, SLL);

theta = (-85:0.1:85)';
sum_beam = zeros(length(theta), 1);
for n = 1:length(theta)
    sv = exp(-1j*2*pi*margin*(0:num_array - 1)'*sind(theta(n)));
    sum_beam(n) = (w_sum.*sv_beam)'*sv;
end

plot(theta, 20*log10(abs(sum_beam)/max(abs(sum_beam))))
grid on
axis([-90, 90, -50, 0])
set(gca, 'XTICK', -90: 15: 90)
set(gca, 'YTICK', -50: 10: 0)

function w_coef = DE_Dif(P, N, w_sum, SLL)
    CR = 0.7;
    F = 0.5;
    f_th = 1e-4;
    MAX_ITR = 1000;
    POP_SIZE = 8;

    w_coef = zeros(P + N, 1);
    Pop = rand(P + N, POP_SIZE);
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
            if newCost < thisCost
                Pop(:, col) = newChromo;
            end
        end
        
    end
    
end

function loss = CostF(w_coef, w_sum, SLL, P)
    coef_sum = w_sum(length(w_sum)/2 + 1:end);
    w_dif = zeros(length(coef_sum), 1);
    for n = P + 1:length(w_coef)
        delta = w_coef(n) == (1:P)';
        w_dif(n - P) = coef_sum(n - P)*sum(delta.*w_coef(1:P));
    end
    w_dif = [-flip(w_dif); w_dif];
    loss = 0;
end