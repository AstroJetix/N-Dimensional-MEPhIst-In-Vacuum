function [E, Accept, M, S] = Ising_Base_Pre(N_x, N_y, J, B, mu_1, T, N_Trial, S_0, E_0, M_0)
% Функция, возвращающая мгновенную энергию системы E [Дж], число принятия
% Accept, мгновенный полный магнитный момент M [А * м^2] и финальную 
% конфигурацию системы
% Функция примерно на 15% быстрее обычной, но имеет ограничение:
% не может работать при слишком большом количестве итераций. 
% N_x - количество частиц по оси X
% N_y - количество частиц по оси Y
% J [Дж] - обменная энергия
% B [Тл] - индукция магнитного поля
% mu_1 [А * м^2] - магнитный момент одной частицы
% T [К] - температура системы
% NTrial - количество испытаний на один набор параметров
% S_0 - начальная матрица конфигурации системы
% E_0 [Дж] - начальная энергия системы
% M_0 [А * м^2] - начальный полный момент системы
k_B = 1.380 * 1e-23; % [Дж / К] Постоянная Больцмана
M = zeros(N_Trial * N_x * N_y + 1, 1);
E = zeros(N_Trial * N_x * N_y + 1, 1);
s = S_0;
E(1) = E_0;
M(1) = M_0;
Accept = 0;
k = 2;
de_mat = [2 * B * mu_1 - 8 * J, 2 * B * mu_1 - 4 * J, 2 * B * mu_1, 2 * B * mu_1 + 4 * J, 2 * B * mu_1 + 8 * J];
exp_mat = exp(-de_mat ./ (k_B * T));
rand_arr = rand(N_Trial * N_x * N_y, 3, "single");
for i=1:N_Trial
    for j=1:N_x*N_y
        I_x = floor(N_x * rand_arr(k - 1, 1) + 1);
        I_y = floor(N_y * rand_arr(k - 1, 2) + 1);
        center = double(s(I_x, I_y));
        up = s(mod(I_x, N_x) + 1, I_y);
        down = s(mod(I_x - 2, N_x) + 1, I_y);
        right = s(I_x, mod(I_y, N_y) + 1);
        left = s(I_x, mod(I_y - 2, N_y) + 1);
        spin_sum = double(up + down + right + left);
        de = center * de_mat(floor(spin_sum / 2 + 3));
        exp_de = exp_mat(floor(spin_sum / 2 + 3)) ^ (center);
        if ((de <= 0) || (rand_arr(k - 1, 3) <= exp_de))
            s(I_x, I_y) = -center;
            Accept = Accept + 1;
            E(k) = E(k - 1) + de;
            M(k) = M(k - 1) - 2 * mu_1 * center;
        else
            E(k) = E(k - 1);
            M(k) = M(k - 1);
        end
        k = k + 1;
    end
end
S = s;
Accept = Accept / (N_Trial * N_x * N_y);