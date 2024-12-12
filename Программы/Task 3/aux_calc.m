syms x y z l q g bx by bz m n t
assume(x, "real")
assume(y, "real")
assume(z, "real")
assume(l, "real")
assume(q, "real")
assume(g, "real")
assume(bx, "real")
assume(by, "real")
assume(bz, "real")
p = sqrt(x^2 + y^2);
r = [x; y];
b = [0.9999; 0]; % Безразмерная скорость
g = 1 / sqrt(1 - norm(b)^2);
q_num = 4.8*1e-10; % [заряд СГС] заряд пробной частицы
lambda_num = 1.5 * 1e-5; % [см] дебаевский радиус
rho_q = -q / (4 * pi) * exp(-g * (p - dot(r, b)) / l) / (l^2 * (p - dot(r, b)));
rho_q = subs(rho_q, q, q_num);
rho_q = subs(rho_q, l, lambda_num);
fsurf(rho_q, [-0.0005 0.0005 -0.0005 0.0005], 'MeshDensity', 100)
%zlim([-100 100])
%colorbar