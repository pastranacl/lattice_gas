dtheta = pi/8.0;
theta = linspace(dtheta, 2*pi - dtheta, 100);

r0 = 1000.0;
A = 0.03;

rsq = r0.*r0.*(1 + A.*sin(6.*theta)).^2;
r = sqrt(rsq);
x = r.*cos(theta);
y = r.*sin(theta);

scatter(real(x), real(y), 10, 'filled')
xlim([-r0*(1+02), r0*(1+02)])
ylim([-r0*(1+02), r0*(1+02)])
