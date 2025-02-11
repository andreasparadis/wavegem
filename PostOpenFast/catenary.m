clear all

d = 100; % Depth
w = 113.35; % Specific weight
b = 0; % Angle initial
z_att = -21.5; % Attachment point
R = 17.5; % Platform radius
X = 250 % Horizontal length of line

Thi(1) = 1e5;

for i = 2:20
%  Thi(i) = (z_att+d)*w / (cosh(w/Thi(i-1)*X)-1);
  Thi(i) = w * X * 1/acosh((z_att+Thi(i-1)/w+d)*w/Thi(i-1));
end

Th = Thi(end) % Tension horizontal
L = Th/w * sinh(w/Th*X)  % Chain length

x = linspace(0,X,1000);
% z = -d + (Th/w * cosh( w/Th.*x + asinh(tan(b))) - Th/w*sqrt(1+tan(b)^2));
z = Th/w * (cosh(w/Th.*x) -1) - d;

ATT = z(end) % Attachment point

figure()
plot(x,z)
axis equal

figure()
plot(1:20, Thi)

A1 = (X+R)*[cos(2*pi/3) sin(2*pi/3)]
A2 = (X+R)*[cos(4*pi/3) sin(4*pi/3)]
A3 = (X+R)*[1 0]

F1 = R*[cos(2*pi/3) sin(2*pi/3)]
F2 = R*[cos(4*pi/3) sin(4*pi/3)]
F3 = R*[1 0]

