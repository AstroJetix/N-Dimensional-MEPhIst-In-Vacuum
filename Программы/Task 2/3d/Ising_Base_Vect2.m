function [E, Accept, M, S] = Ising_Base_Vect2(N_x, N_y, N_z, J, B, mu_1, T, N_Trial, ...
                                        S_0, E_0, M_0)
% Функция, возвращающая мгновенную энергию системы E [Дж], число принятия
% Accept, мгновенный полный магнитный момент M [А * м^2] и финальную 
% конфигурацию системы
% N_x - количество частиц по оси X
% N_y - количество частиц по оси Y
% N_z - количество частиц по оси Z
% J [Дж] - обменная энергия
% B [Тл] - индукция магнитного поля
% mu_1 [А * м^2] - магнитный момент одной частицы
% T [К] - температура системы
% NTrial - количество испытаний на один набор параметров
% S_0 - начальная матрица конфигурации системы
% E_0 [Дж] - начальная энергия системы
% M_0 [А * м^2] - начальный полный момент системы
k_B = 1.380 * 1e-23; % [Дж / К] Постоянная Больцмана
M = zeros(N_Trial * N_x * N_y * N_z, 1);
E = zeros(N_Trial * N_x * N_y * N_z, 1);
E(1) = E_0;
M(1) = M_0;
Accept = 0;
k = 2;
% Создание циклической матрицы
Circ_s = zeros(N_x + 2, N_y + 2, N_z + 2);
Circ_s(2:N_x+1, 2:N_y+1, 2:N_z+1) = S_0;
Circ_s(1, 2:N_y+1, 2:N_z+1) = S_0(N_x, :, :);
Circ_s(N_x+2, 2:N_y+1, 2:N_z+1) = S_0(1, :, :);
Circ_s(2:N_x+1, 1, 2:N_z+1) = S_0(:, N_y, :);
Circ_s(2:N_x+1, N_y+2, 2:N_z+1) = S_0(:, 1, :);
Circ_s(2:N_x+1, 2:N_y+1, 1) = S_0(:, :, N_z);
Circ_s(2:N_x+1, 2:N_y+1, N_z+2) = S_0(:, :, 1);
Circ_s = int8(Circ_s);
% Создание маски для быстрого вычисления суммы спинов
mask = int8(zeros(N_x + 2, N_y + 2, N_z + 2));
for i=1:N_Trial
    for j=1:N_x*N_y*N_z
        I_x = floor(N_x * rand() + 1) + 1;
        I_y = floor(N_y * rand() + 1) + 1;
        I_z = floor(N_z * rand() + 1) + 1;
        center = Circ_s(I_x, I_y, I_z);
        % Наполнение маски
        if (N_x ~= 1)
            mask(I_x-1, I_y, I_z) = 1;
            mask(I_x+1, I_y, I_z) = 1;
        end
        if (N_y ~= 1)
            mask(I_x, I_y-1, I_z) = 1;
            mask(I_x, I_y+1, I_z) = 1;
        end
        if (N_z ~= 1)
            mask(I_x, I_y, I_z-1) = 1;
            mask(I_x, I_y, I_z+1) = 1;
        end
        spin_sum = sum(Circ_s.*mask, "all");
        % Обнуляем маску
        mask(I_x-1, I_y, I_z) = 0;
        mask(I_x+1, I_y, I_z) = 0;
        mask(I_x, I_y-1, I_z) = 0;
        mask(I_x, I_y+1, I_z) = 0;
        mask(I_x, I_y, I_z-1) = 0;
        mask(I_x, I_y, I_z+1) = 0;
        de = 2 * double(center) * B * mu_1 + 2 * J * double(center * spin_sum); 
        if ((de <= 0) || (rand() <= exp(-de ./ (k_B * T))))
            Circ_s(I_x, I_y, I_z) = -center;
            % Восстановление цикличности матрицы
            if ((N_x ~= 1) && (I_x == 2))
                Circ_s(N_x + 2, I_y, I_z) = -center;
            end
            if ((N_x ~= 1) && (I_x == N_x+1))
                Circ_s(1, I_y, I_z) = -center;
            end
            if ((N_y ~= 1) && (I_y == 2))
                Circ_s(I_x, N_y + 2, I_z) = -center;
            end
            if ((N_y ~= 1) && (I_y == N_y+1))
                Circ_s(I_x, 1, I_z) = -center;
            end
            if ((N_z ~= 1) && (I_z == 2))
                Circ_s(I_x, I_y, N_z+2) = -center;
            end
            if ((N_z ~= 1) && (I_z == N_z+1))
                Circ_s(I_x, I_y, 1) = -center;
            end
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
S = Circ_s(2:N_x+1, 2:N_y+1, 2:N_z+1);
Accept = Accept / (N_Trial * N_x * N_y * N_z);