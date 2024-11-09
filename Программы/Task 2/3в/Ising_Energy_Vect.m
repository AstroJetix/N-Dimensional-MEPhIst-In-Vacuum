function z = Ising_Energy_Vect(N_x, N_y, N_z, S, J, B, mu_1)
% Функция, возвращающая энергию данной конфигурации системы
% N_x - количество частиц по оси X
% N_y - количество частиц по оси Y
% N_z - количество частиц по оси Z
% S - матрица конфигурации системы
% J [Дж] - обменная энергия
% B [Тл] - индукция магнитного поля
% mu_1 [А * м^2] - магнитный момент одной частицы
E = 0;
% Составление циклической матрицы
%S_Circ = int8(zeros(N_x + 1, N_y + 1, N_z + 1));
S_Circ = zeros(N_x + 1, N_y + 1, N_z + 1);
S_Circ(1:N_x, 1:N_y, 1:N_z) = S;
S_Circ(N_x+1, 1:N_y, 1:N_z) = S(1, 1:N_y, 1:N_z);
S_Circ(1:N_x, N_y+1, 1:N_z) = S(1:N_x, 1, 1:N_z);
S_Circ(1:N_x, 1:N_y, N_z+1) = S(1:N_x, 1:N_y, 1);
% Подсчет энергий по каждой оси 
if (N_x ~= 1)
    E = E - J * sum(S_Circ(2:N_x+1, 1:N_y, 1:N_z).*S, "all");
end
if (N_y ~= 1)
    E = E - J * sum(S_Circ(1:N_x, 2:N_y+1, 1:N_z).*S, "all");
end
if (N_z ~= 1)
    E = E - J * sum(S_Circ(1:N_x, 1:N_y, 2:N_z+1).*S, "all");
end
E = E - mu_1 * B * sum(S, "all");
z = E;