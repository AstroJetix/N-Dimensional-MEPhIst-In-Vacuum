beta = [0.5; 0];
g = 1 / sqrt(1 - norm(beta)^2)
q = 1;
lambda = 1;
dx = 0.1;
xmin = -5 + dx / 2;
xmax = 5 - dx / 2;
[X, Y] = meshgrid(xmin:dx:xmax, xmin:dx:xmax);
r = sqrt(X.^2 + Y.^2);
field_x = q * (X - beta(1) * r) ./ (g * (r - beta(1) * X - beta(2) * Y).^2)...
        .* exp(-g * (r - beta(1) * X - beta(2) * Y) / lambda) .* (1 / lambda + 1 ./ (g * (r - beta(1) * X - beta(2) * Y)));
field_y = q * (Y - beta(2) * r) ./ (g * (r - beta(1) * X - beta(2) * Y).^2)...
        .* exp(-g * (r - beta(1) * X - beta(2) * Y) / lambda) .* (1 / lambda + 1 ./ (g * (r - beta(1) * X - beta(2) * Y)));

quiver(X, Y, field_x, field_y)
D = divergence(X, Y, field_x, field_y);
%D(D > 0) = 0;
hold on
contourf(X, Y, D)
colorbar
hold off