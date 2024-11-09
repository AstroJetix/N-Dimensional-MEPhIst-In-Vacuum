function [E, Accept, M, S] = Ising_Base(N_x, N_y, N_z, J, B, mu_1, T, N_Trial, ...
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
s = S_0;
E(1) = E_0;
M(1) = M_0;
Accept = 0;
k = 2;
for i=1:N_Trial
    for j=1:N_x*N_y*N_z
        I_x = floor(N_x * rand() + 1);
        I_y = floor(N_y * rand() + 1);
        I_z = floor(N_z * rand() + 1);
        center = s(I_x, I_y, I_z);
        fwd = mod(I_x, N_x) + 1;
        bck = mod(I_x - 2, N_x) + 1;
        right = mod(I_y, N_y) + 1;
        left = mod(I_y - 2, N_y) + 1;
        up = mod(I_z, N_z) + 1;
        down = mod(I_z - 2, N_z) + 1;
        de = 2 * double(s(I_x, I_y, I_z)) * B * mu_1; 
        if (N_x ~= 1)
            %de = de + 2 * J * double(s(I_x, I_y, I_z) * (s(fwd, I_y, I_z) + s(bck, I_y, I_z)));
            fwd = s(fwd, I_y, I_z);
            bck = s(bck, I_y, I_z);
            de = de + 2 * J * center * (fwd + bck);
        end
        if (N_y ~= 1)
            %de = de + 2 * J * double(s(I_x, I_y, I_z) * (s(I_x, right, I_z) + s(I_x, left, I_z)));
            left = s(I_x, left, I_z);
            right = s(I_x, right, I_z);
            de = de + 2 * J * center * (left + right);
        end
        if (N_z ~= 1)
            %de = de + 2 * J * double(s(I_x, I_y, I_z) * (s(I_x, I_y, up) + s(I_x, I_y, down)));
            de = de + 2 * J * s(I_x, I_y, I_z) * (s(I_x, I_y, up) + s(I_x, I_y, down));
            up = s(I_x, I_z, up);
            down = s(I_x, I_y, down);
            de = de + 2 * J * center * (up + down);
        end
        if ((de <= 0) || (rand() <= exp(-de ./ (k_B * T))))
            s(I_x, I_y, I_z) = -s(I_x, I_y, I_z);
            Accept = Accept + 1;
            E(k) = E(k - 1) + de;
            M(k) = mu_1 * sum(s, "all");
        else
            E(k) = E(k - 1);
            M(k) = M(k - 1);
        end
        k = k + 1;
    end
end
S = s;
Accept = Accept / (N_Trial * N_x * N_y * N_z);