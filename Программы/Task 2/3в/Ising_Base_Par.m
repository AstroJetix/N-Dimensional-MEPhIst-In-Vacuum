function [E, Accept, M, S] = Ising_Base_Par(N_x, N_y, N_z, J, B, mu_1, T, N_Trial, ...
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
M_int = zeros(N_x * N_y * N_z, 1);
E_int = zeros(N_x * N_y * N_z, 1);
s = S_0;
E(1) = E_0;
M(1) = M_0;
Accept = 0;
Circ_s = zeros(N_x + 2, N_y + 2, N_z + 2);
Circ_s(2:N_x+1, 2:N_y+1, 2:N_z+1) = S_0;
Circ_s(1, 2:N_y+1, 2:N_z+1) = S_0(N_x, :, :);
Circ_s(N_x+2, 2:N_y+1, 2:N_z+1) = S_0(1, :, :);
Circ_s(2:N_x+1, 1, 2:N_z+1) = S_0(:, N_y, :);
Circ_s(2:N_x+1, N_y+2, 2:N_z+1) = S_0(:, 1, :);
Circ_s(2:N_x+1, 2:N_y+1, 1) = S_0(:, :, N_z);
Circ_s(2:N_x+1, 2:N_y+1, N_z+2) = S_0(:, :, 1);
mask = zeros(3, 3, 3);
if (N_x ~= 1)
    mask(1, 2, 2) = 1;
    mask(3, 2, 2) = 1;
end
if (N_y ~= 1)
    mask(2, 1, 2) = 1;
    mask(2, 3, 2) = 1;
end
if (N_z ~= 1)
    mask(2, 2, 1) = 1;
    mask(2, 2, 3) = 1;
end
for i=1:N_Trial
    parfor j=1:N_x*N_y*N_z
        I_x = floor(N_x * rand() + 1);
        I_y = floor(N_y * rand() + 1);
        I_z = floor(N_z * rand() + 1);
        
        de = 2 * double(cent) * B * mu_1; 
        
        if ((de <= 0) || (rand() <= exp(-de ./ (k_B * T))))
            s(I_x, I_y, I_z) = -cent;
            Accept = Accept + 1;
        end
        E_int(j) = Ising_Energy_Vect(N_x, N_y, N_z, s, J, B, mu_1);
        M_int(j) = mu_1 * sum(s, "all");
    end
end
S = s;
Accept = Accept / (N_Trial * N_x * N_y * N_z);