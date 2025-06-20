function [xK, yK, zK, xdotK, ydotK, zdotK] = computeGeneralSolutionKechichian(J2, J3, R, a, e, INC, Omega, omega, theta_0, X0, n, t)

% Vector inicial de estados
x0  = X0(1);
y0  = X0(2);
z0  = X0(3);
dx0 = X0(4);
dy0 = X0(5);
dz0 = X0(6);

% Constantes
nt = n*t;
cnt = cos(nt);
snt = sin(nt);

si = sin(INC);
ci = cos(INC);
sw = sin(omega);
cw = cos(omega);

delta = 2*(theta_0 - e*(1 + sqrt(1 - e^2))*sin(theta_0));
sd2 = sin(delta/2);
cd2 = sin(delta/2);
theta_c = omega + theta_0 - 2*e*sd2;
AOL = theta_c + nt + 2*e*(sd2*cnt + cd2*snt);

f0    =   - e.*cos(theta_c - delta./2);
f(1)  =     - sin(theta_c);
f(2)  =       cos(theta_c);
f(3)  =   - e.*sin(theta_c + delta./2);
f(4)  =     e.*cos(theta_c + delta./2);
f(5)  =   2.*e.*sin(2.*theta_c - delta./2);
f(6)  = - 2.*e.*cos(2.*theta_c - delta./2);
f(7)  =     - sin(2.*theta_c);
fp(7) = - 2.*e.*sin(2.*theta_c + delta./2);
f(8)  =       cos(2.*theta_c);
fp(8) =   2.*e.*cos(2.*theta_c + delta./2);

g0    =   - e.*sin(theta_c - delta./2);
g(1)  =       cos(theta_c);
g(2)  =       sin(theta_c);
g(3)  =     e.*cos(theta_c + delta./2);
g(4)  =     e.*sin(theta_c + delta./2);
g(5)  = - 2.*e.*cos(2.*theta_c - delta./2);
g(6)  = - 2.*e.*sin(2.*theta_c - delta./2);
g(7)  =       cos(2.*theta_c);
gp(7) =   2.*e.*cos(2.*theta_c + delta./2);
g(8)  =       sin(2.*theta_c);
gp(8) =   2.*e.*sin(2.*theta_c + delta./2);

ac(1)  = - 3.*n.^2.*J2.*R.^2.*si.^2./(a.*(1 - e.^2).^4);
ap(1)  = - 3.*n.^2.*J2.*R.^2./(2.*a.*(1 - e.^2).^4);
app(1) = - 3.*n.^2.*J2.*R.^2.*si.*ci./(1 - e.^2).^4;
ac(2)  =   cw.*f0 +   sw.*g0;
ac(3)  = cw.*f(1) + sw.*g(1);
ac(4)  = cw.*f(2) + sw.*g(2);
ac(5)  = cw.*f(3) + sw.*g(3);
ac(6)  = cw.*f(4) + sw.*g(4);

b0    = ac(1).*e.*(ac(3).*g(5) + ac(4).*g(6) + ac(5).*g(7) + ac(6).*g(8));
b(1)  = ac(1).*( g(5)./2 + 2.*e.*ac(2).*g(5) -   e.*ac(6).*g(5) +   e.*ac(5).*g(6) +  e.*ac(4).*g(7)  -   e.*ac(3).* g(8) + e.*ac(6).*gp(7) - e.*ac(5).*gp(8));
bp(1) = ac(1).*( g(6)./2 +   e.*ac(5).*g(5) + 2.*e.*ac(2).*g(6) +   e.*ac(6).*g(6) +  e.*ac(3).*g(7)  +   e.*ac(4).* g(8) + e.*ac(5).*gp(7) + e.*ac(6).*gp(8));
b(2)  = ac(1).*( g(7)./2 +   e.*ac(4).*g(5) +   e.*ac(3).*g(6) + 2.*e.*ac(2).*g(7) + e.*ac(4).*gp(7) -   e.*ac(3).*gp(8));
bp(2) = ac(1).*( g(8)./2 -   e.*ac(3).*g(5) +   e.*ac(4).*g(6) + 2.*e.*ac(2).*g(8) + e.*ac(3).*gp(7) +   e.*ac(4).*gp(8));
b(3)  = ac(1).*(gp(7)./2 +   e.*ac(6).*g(5) +   e.*ac(5).*g(6) +   e.*ac(4).*g(7) +  e.*ac(3).*g(8) + 2.*e.*ac(2).*gp(7));
bp(3) = ac(1).*(gp(8)./2 -   e.*ac(5).*g(5) +   e.*ac(6).*g(6) -   e.*ac(3).*g(7) +  e.*ac(4).*g(8) + 2.*e.*ac(2).*gp(8));
b(4)  = ac(1).*e.*( ac(6).*g(7) + ac(5).*g(8) + ac(4).*gp(7) + ac(3).*gp(8));
bp(4) = ac(1).*e.*(-ac(5).*g(7) + ac(6).*g(8) - ac(3).*gp(7) + ac(4).*gp(8));
b(5)  = ac(1).*e.*(ac(6).*gp(7) + ac(5).*g(8));
bp(5) = ac(1).*e.*(ac(6).*gp(8) - ac(5).*g(7));

c0    = ap(1).*(1 - 3./2.*si.^2  + ac(2).*e.*(4 - 6.*si.^2) + 6.*e.*si.^2.*( ac(3)./2.*f(5) + ac(4)./2.*f(6) + ac(5)./2.*f(7) + ac(6)./2.*f(8)));
c(1)  = ap(1).*(3./2.*si.^2.*f(5) + ac(3).*e.*(4 - 6.*si.^2) + 6.*e.*si.^2.*( ac(2)  .*f(5) - ac(6)./2.*f(5) + ac(5)./2.*f(6) + ac(4)./2.*f(7) - ac(3)./2.*f(8) + ac(6)./2.*fp(7) - ac(5)./2.*fp(8)));
cp(1) = ap(1).*(3./2.*si.^2.*f(6) + ac(4).*e.*(4 - 6.*si.^2) + 6.*e.*si.^2.*( ac(5)./2.*f(5) + ac(2)  .*f(6) + ac(6)./2.*f(6) + ac(3)./2.*f(7) + ac(4)./2.*f(8) + ac(5)./2.*fp(7) + ac(6)./2.*fp(8)));
c(2)  = ap(1).*(3./2.*si.^2.*f(7) + ac(5).*e.*(4 - 6.*si.^2) + 6.*e.*si.^2.*( ac(4)./2.*f(5) + ac(3)./2.*f(6) + ac(2)  .*f(7) + ac(4)./2.*fp(7) - ac(3)./2.*fp(8)));
cp(2) = ap(1).*(3./2.*si.^2.*f(8) + ac(6).*e.*(4 - 6.*si.^2) + 6.*e.*si.^2.*(-ac(3)./2.*f(5) + ac(4)./2.*f(6) + ac(2)  .*f(8) + ac(3)./2.*fp(7) + ac(4)./2.*fp(8)));
c(3)  = 3./2.*si.^2.*ap(1).*fp(7) + ap(1).*6.*e.*si.^2.*(  ac(6)./2.*f(5) + ac(5)./2.*f(6) + ac(4)./2.*f(7) + ac(3)./2.*f(8) + ac(2).*fp(7));
cp(3) = 3./2.*si.^2.*ap(1).*fp(8) + ap(1).*6.*e.*si.^2.*(- ac(5)./2.*f(5) + ac(6)./2.*f(6) - ac(3)./2.*f(7) + ac(4)./2.*f(8) + ac(2).*fp(8));
c(4)  = ap(1).*6.*e.*si.^2.*(  ac(6)./2.*f(7) + ac(5)./2.*f(8) + ac(4)./2.*fp(7) + ac(3)./2.*fp(8));
cp(4) = ap(1).*6.*e.*si.^2.*(- ac(5)./2.*f(7) + ac(6)./2.*f(8) - ac(3)./2.*fp(7) + ac(4)./2.*fp(8));
c(5)  = 3.*ap(1).*e.*si.^2.*(ac(6).*fp(7) + ac(5).*fp(8));
cp(5) = 3.*ap(1).*e.*si.^2.*(ac(6).*fp(8) - ac(5).*fp(7));

d0    = app(1).*(g0   + 4.*e.*ac(2).*g0   + 2.*e.*ac(3).*g(1) + 2.*e.*ac(4).*g(2) + 2.*e.*ac(5).*g(3) + 2.*e.*ac(6).*g(4));
d(1)  = app(1).*(g(1) + 4.*e.*ac(2).*g(1) + 4.*e.*ac(3).*g0   - 2.*e.*ac(3).*g(4) + 2.*e.*ac(4).*g(3) + 2.*e.*ac(5).*g(2) - 2.*e.*ac(6).*g(1));
dp(1) = app(1).*(g(2) + 4.*e.*ac(2).*g(2) + 2.*e.*ac(3).*g(3) + 2.*e.*ac(4).*g(4) + 4.*e.*ac(4).*g0   + 2.*e.*ac(5).*g(1) + 2.*e.*ac(6).*g(2));
d(2)  = app(1).*(g(3) + 4.*e.*ac(2).*g(3) + 2.*e.*ac(3).*g(2) + 2.*e.*ac(4).*g(1) + 4.*e.*ac(5).*g0);
dp(2) = app(1).*(g(4) + 4.*e.*ac(2).*g(4) - 2.*e.*ac(3).*g(1) + 2.*e.*ac(4).*g(2) + 4.*e.*ac(6).*g0);
d(3)  = app(1).*2.*e.*(  ac(3).*g(4) + ac(4).*g(3) + ac(5).*g(2) + ac(6).*g(1));
dp(3) = app(1).*2.*e.*(- ac(3).*g(3) + ac(4).*g(4) - ac(5).*g(1) + ac(6).*g(2));
d(4)  = app(1).*2.*e.*(  ac(5).*g(4) + ac(6).*g(3));
dp(4) = app(1).*2.*e.*(- ac(5).*g(3) + ac(6).*g(4));

ep1 = dz0./n;
ep2 = z0;
kp1 = dx0./n;
kp2 = x0;
K   = dy0 + 2.*n.*x0 + sum(b(1:5) ./ ((1:5)*n));
Kp  = y0 - 2.*kp1 + 4.*b0./n.^2 ...
      - sum(2./((2:5)*n^2) .* (2.*bp(2:5)./(2:5) + c(2:5)))  ...
      + sum(bp(1:5) ./ ((1:5).^2 * n^2)) ...
      - 2./n.^2.*(2.*bp(1) + c(1));

f(9)  =   3.*e.*sin(3.*theta_c - delta./2);
f(10) = - 3.*e.*cos(3.*theta_c - delta./2);
f(11) =     - sin(3.*theta_c);
f(12) =       cos(3.*theta_c);
f(13) = - 3.*e.*sin(3.*theta_c + delta./2);
f(14) =   3.*e.*cos(3.*theta_c + delta./2);

g(9)  = - 3.*e.*cos(3.*theta_c - delta./2);
g(10) = - 3.*e.*sin(3.*theta_c - delta./2);
g(11) =       cos(3.*theta_c);
g(12) =       sin(3.*theta_c);
g(13) =   3.*e.*cos(3.*theta_c + delta./2);
g(14) =   3.*e.*sin(3.*theta_c + delta./2);

I(1) = - n.^2.*J3.*R.^3./(2.*a.^2.*(1 - e.^2).^5);
I(2) =   15./4.*si.^3 - 3.*si;
I(3) = - 15./4.*si.^3;

m0    = I(1).*I(2).*f0   + 5.*I(1).*I(2).*e.*(ac(2).*f0 + 1./2.*ac(3).*f(1) + ac(4)./2.*f(2) + ac(5)./2.*f(3) + ac(6)./2.*f(4)) ...
        + 5.*I(1).*I(3).*e.*(ac(5)./2.*f(9) + ac(6)./2.*f(10));
m(1)  = I(1).*I(2).*f(1) + 5.*I(1).*I(2).*e.*(ac(2).*f(1) + ac(3).*f0 - ac(3)./2.*f(4) + ac(4)./2.*f(3) + ac(5)./2.*f(2) - ac(6)./2.*f(1)) ...
        + 5.*I(1).*I(3).*e.*(- ac(3)./2.*f(10) + ac(4)./2.*f(9) - ac(5)./2.*f(12) + ac(6)./2.*f(11));
mp(1) = I(1).*I(2).*f(2) + 5.*I(1).*I(2).*e.*(ac(2).*f(2) + ac(3)./2.*f(3) + ac(4).*f0 + ac(4)./2.*f(4) + ac(5)./2.*f(1) + ac(6)./2.*f(2)) ...
        + 5.*I(1).*I(3).*e.*(ac(3)./2.*f(9) + ac(4)./2.*f(10) + ac(5)./2.*f(11) + ac(6)./2.*f(12));
m(2)  = I(1).*I(2).*f(3) + I(1).*I(3).*f(9) + 5.*I(1).*I(2).*e.*(ac(2).*f(3) + ac(3)./2.*f(2) + ac(4)./2.*f(1) + ac(5).*f0) ...
        + 5.*I(1).*I(3).*e.*(ac(2).*f(9) - ac(3)./2.*f(12) + ac(4)./2.*f(11) - ac(5)./2.*f(14) + ac(6)./2.*f(13));
mp(2) = I(1).*I(2).*f(4) + I(1).*I(3).*f(10) + 5.*I(1).*I(2).*e.*(ac(2).*f(4) - ac(3)./2.*f(1) + ac(4)./2.*f(2) + ac(6).*f0) ...
        + 5.*I(1).*I(3).*e.*(ac(2).*f(10) + ac(3)./2.*f(11) + ac(4)./2.*f(12) + ac(5)./2.*f(13) + ac(6)./2.*f(14));
m(3)  = I(1).*I(3).*f(11) + 5.*I(1).*I(2).*e.*(ac(3)./2.*f(4) + ac(4)./2.*f(3) + ac(5)./2.*f(2) + ac(6)./2.*f(1)) ...
        + 5.*I(1).*I(3).*e.*(ac(2).*f(11) + ac(3)./2.*f(10) - ac(3)./2.*f(14) + ac(4)./2.*f(9) + ac(4)./2.*f(13));
mp(3) = I(1).*I(3).*f(12) + 5.*I(1).*I(2).*e.*(-ac(3)./2.*f(3) + ac(4)./2.*f(4) - ac(5)./2.*f(1) + ac(6)./2.*f(2)) ...
        + 5.*I(1).*I(3).*e.*(ac(2).*f(12) - ac(3)./2.*f(9) + ac(3)./2.*f(13) + ac(4)./2.*f(10) + ac(4)./2.*f(14));
m(4)  = I(1).*I(3).*f(13) + 5.*I(1).*I(2).*e.*(ac(5)./2.*f(4) + ac(6)./2.*f(3)) ...
        + 5.*I(1).*I(3).*e.*(ac(2).*f(13) + ac(3)./2.*f(12) + ac(4)./2.*f(11) + ac(5)./2.*f(10) + ac(6)./2.*f(9));
mp(4) = I(1).*I(3).*f(14) + 5.*I(1).*I(2).*e.*(-ac(5)./2.*f(3) + ac(6)./2.*f(4)) ...
        + 5.*I(1).*I(3).*e.*(ac(2).*f(14) - ac(3)./2.*f(11) + ac(4)./2.*f(12) - ac(5).*f(9) + ac(6)./2.*f(10));
m(5)  = 5.*I(1).*I(3).*e.*( ac(3)./2.*f(14) + ac(4)./2.*f(13) + ac(5)./2.*f(12) + ac(6)./2.*f(11));
mp(5) = 5.*I(1).*I(3).*e.*(-ac(3)./2.*f(13) + ac(4)./2.*f(14) - ac(5)./2.*f(11) + ac(6)./2.*f(12));
m(6)  = 5.*I(1).*I(3).*e.*( ac(5)./2.*f(14) + ac(6)./2.*f(13));
mp(6) = 5.*I(1).*I(3).*e.*(-ac(5)./2.*f(13) + ac(6)./2.*f(14));

p(1) = - n.^2.*J3.*R.^3./(a.^2.*(1 - e.^2).^5);
p(2) =   6.*si - 15./2.*si.^3;
p(3) =   5./2.*si.^3;

r0    = p(1).*p(2).*g0 + 5.*p(1).*p(2).*e.*(ac(2).*g0 + ac(3)./2.*g(1) + ac(4)./2.*g(2) + ac(5)./2.*g(3) + ac(6)./2.*g(4)) ...
        + 5.*p(1).*p(3).*e.*(ac(5)./2.*g(9) + ac(6)./2.*g(10));
r(1)  = p(1).*p(2).*g(1) + 5.*p(1).*p(2).*e.*(ac(2).*g(1) + ac(3).*g0 - ac(3)./2.*g(4) + ac(4)./2.*g(3) + ac(5)./2.*g(2) - ac(6)./2.*g(1)) ...
        + 5.*p(1).*p(3).*e.*(-ac(3)./2.*g(10) + ac(4)./2.*g(9) - ac(5)./2.*g(12) + ac(6)./2.*g(11));
rp(1) = p(1).*p(2).*g(2) + 5.*p(1).*p(2).*e.*(ac(2).*g(2) + ac(3)./2.*g(3) + ac(4).*g0 + ac(4)./2.*g(4) + ac(5)./2.*g(1) + ac(6)./2.*g(2)) ...
        + 5.*p(1).*p(3).*e.*(ac(3)./2.*g(9) + ac(4)./2.*g(10) + ac(5)./2.*g(11) + ac(6)./2.*g(12));
r(2)  = p(1).*p(2).*g(3) + p(1).*p(3).*g(9) + 5.*p(1).*p(2).*e.*(ac(2).*g(3) + ac(3)./2.*g(2) + ac(4)./2.*g(1) + ac(5).*g0) ...
        + 5.*p(1).*p(3).*e.*(ac(2).*g(9) - ac(3)./2.*g(12) + ac(4)./2.*g(11) - ac(5)./2.*g(14) + ac(6)./2.*g(13));
rp(2) = p(1).*p(2).*g(4) + p(1).*p(3).*g(10) + 5.*p(1).*p(2).*e.*(ac(2).*g(4) - ac(3)./2.*g(1) + ac(4)./2.*g(2) + ac(6).*g0) ...
        + 5.*p(1).*p(3).*e.*(ac(2).*g(10) + ac(3)./2.*g(11) + ac(4)./2.*g(12) + ac(5)./2.*g(13) + ac(6)./2.*g(14));
r(3)  = p(1).*p(3).*g(11) + 5.*p(1).*p(2).*e.*(ac(3)./2.*g(4) + ac(4)./2.*g(3) + ac(5)./2.*g(2) + ac(6)./2.*g(1)) ...
        + 5.*p(1).*p(3).*e.*(ac(2).*g(11) + ac(3)./2.*g(10) - ac(3)./2.*g(14) + ac(4)./2.*g(9) + ac(4)./2.*g(13));
rp(3) = p(1).*p(3).*g(12) + 5.*p(1).*p(2).*e.*(-ac(3)./2.*g(3) + ac(4)./2.*g(4) - ac(5)./2.*g(1) + ac(6)./2.*g(2)) ...
        + 5.*p(1).*p(3).*e.*(ac(2).*g(12) - ac(3)./2.*g(9) + ac(3)./2.*g(13) + ac(4)./2.*g(10) + ac(4)./2.*g(14));
r(4)  = p(1).*p(3).*g(13) + 5.*p(1).*p(2).*e.*(ac(5)./2.*g(4) + ac(6)./2.*g(3)) ...
        + 5.*p(1).*p(3).*e.*(ac(2).*g(13) + ac(3)./2.*g(12) + ac(4)./2.*g(11) + ac(5)./2.*g(10) + ac(6)./2.*g(9));
rp(4) = p(1).*p(3).*g(14) + 5.*p(1).*p(2).*e.*(-ac(5)./2.*g(3) + ac(6)./2.*g(4)) ...
        + 5.*p(1).*p(3).*e.*(ac(2).*g(14) - ac(3)./2.*g(11) + ac(4)./2.*g(12) - ac(5)./2.*g(9) + ac(6)./2.*g(10));
r(5)  = 5.*p(1).*p(3).*e.*( ac(3)./2.*g(14) + ac(4)./2.*g(13) + ac(5)./2.*g(12) + ac(6)./2.*g(11));
rp(5) = 5.*p(1).*p(3).*e.*(-ac(3)./2.*g(13) + ac(4)./2.*g(13) - ac(5)./2.*g(11) + ac(6)./2.*g(12));
r(6)  = 5.*p(1).*p(3).*e.*( ac(5)./2.*g(14) + ac(6)./2.*g(13));
rp(6) = 5.*p(1).*p(3).*e.*(-ac(5)./2.*g(13) + ac(6)./2.*g(14));

tc(1) = - n.^2.*J3.*R.^3./(2.*a.^2.*(1 - e.^2).^5);
tc(2) =   15./2.*si.^2.*ci - 3.*ci;
tc(3) = - 15./2.*si.^2.*ci;

u0    = tc(1).*tc(2) + 5.*tc(1).*tc(2).*e.*ac(2) + 5.*tc(1).*tc(3).*e.*(ac(3)./2.*f(5) + ac(4)./2.*f(6) + ac(5)./2.*f(7) + ac(6)./2.*f(8));
u(1)  = tc(1).*tc(3).*f(5) + 5.*tc(1).*tc(2).*e.*ac(3) + ...
        5.*tc(1).*tc(3).*e.*(ac(2).*f(5) - ac(3)./2.*f(8) + ac(4)./2.*f(7) + ac(5)./2.*f(6) - ac(6)./2.*f(5) - ac(3)./2.*fp(8) - ac(5)./2.*fp(8) + ac(6)./2.*fp(7));
up(1) = tc(1).*tc(3).*f(6) + 5.*tc(1).*tc(2).*e.*ac(4) + ...
        5.*tc(1).*tc(3).*e.*(ac(2).*f(6) + ac(3)./2.*f(7) + ac(4)./2.*f(8) + ac(5)./2.*f(5) + ac(6)./2.*f(6) + ac(6)./2.*fp(8) + ac(5)./2.*fp(7) + ac(3)./2.*fp(7));
u(2)  = tc(1).*tc(3).*f(7) + 5.*tc(1).*tc(2).*e.*ac(5) + 5.*tc(1).*tc(3).*e.*(ac(2).*f(7) + ac(3)./2.*f(6) + ac(4)./2.*f(5) + ac(4)./2.*fp(7));
up(2) = tc(1).*tc(3).*f(8) + 5.*tc(1).*tc(2).*e.*ac(6) + 5.*tc(1).*tc(3).*e.*(ac(2).*f(8) - ac(3)./2.*f(5) + ac(4)./2.*f(6) + ac(4)./2.*fp(8));
u(3)  = 5.*tc(1).*tc(3).*e.*(  ac(3)./2.*f(8) + ac(4)./2.*f(7) + ac(5)./2.*f(6) + ac(6)./2.*f(5) + ac(2).*fp(7) + ac(3)./2.*fp(8));
up(3) = 5.*tc(1).*tc(3).*e.*(- ac(3)./2.*f(7) + ac(4)./2.*f(8) - ac(5)./2.*f(5) + ac(6)./2.*f(6) + ac(2).*fp(8) + ac(3)./2.*fp(7));
u(4)  = 5.*tc(1).*tc(3).*e.*(  ac(5)./2.*f(8) + ac(6)./2.*f(7) + ac(4)./2.*fp(7));
up(4) = 5.*tc(1).*tc(3).*e.*(- ac(5)./2.*f(7) + ac(6)./2.*f(8) + ac(4)./2.*fp(8));
u(5)  = 5.*tc(1).*tc(3).*e.*(ac(6)./2.*fp(7) + ac(5)./2.*fp(8));
up(5) = 5.*tc(1).*tc(3).*e.*(ac(6)./2.*fp(8) - ac(5)./2.*fp(7));

K1  = dy0 + 2.*n.*x0 + sum(m(1:6) ./ ((1:6)*n));
Kp1 = y0 - 2.*dx0./n + 4.*m0./n.^2 ...
      - sum(2./((2:6)*n^2) .* (2.*mp(2:6)./(2:6) + r(2:6))) ...
      + sum(mp(1:6) ./ ((1:6).^2 * n^2)) ...
      - 2.*(2.*mp(1) + r(1))./n.^2;

% Solución analítica
% Inicialización de sumas
 sum1 = 0;  sum2 = 0;  sum3 = 0;  sum4 = 0;  sum5 = 0;  
 sum6 = 0;  sum7 = 0;  sum8 = 0;  sum9 = 0; sum10 = 0; 
sum11 = 0; sum12 = 0; sum13 = 0; sum14 = 0; sum15 = 0; 
sum16 = 0; sum17 = 0; sum18 = 0; sum19 = 0; sum20 = 0; 
sum21 = 0; sum22 = 0; sum23 = 0; sum24 = 0; sum25 = 0; 
sum26 = 0; sum27 = 0; sum28 = 0; sum29 = 0; sum30 = 0;
sum31 = 0; sum32 = 0;

for j = 1:5
    sum5  = sum5 + (b(j)/(j^2*n^2)) * sin(j*nt);
    sum6  = sum6 + (bp(j)/(j^2*n^2)) * cos(j*nt);
    sum13 = sum13 + b(j)/(j*n) * cos(j*nt);
    sum14 = sum14 + bp(j)/(j*n) * sin(j*nt);
end
for j = 1:6
    sum21 = sum21 + m(j)/(j^2*n^2) * sin(j*nt);
    sum22 = sum22 + mp(j)/(j^2*n^2) * cos(j*nt);
    sum29 = sum29 + m(j)/(j*n) * cos(j*nt);
    sum30 = sum30 + mp(j)/(j*n) * sin(j*nt);
end
for j = 2:4
    sum7  = sum7 + d(j)*(j*snt - sin(j*nt))/(n^2*(j^2 - 1));
    sum8  = sum8 + dp(j)*(cnt - cos(j*nt))/(n^2*(j^2 - 1));
    sum15 = sum15 + d(j) * (j*cnt - j*cos(j*nt)) / (n*(j^2 - 1));
    sum16 = sum16 + dp(j) * (-snt + j*sin(j*nt)) / (n*(j^2 - 1));
end
for j = 2:5
    sum1  = sum1 + (2*bp(j)/j + c(j)) * (j*snt - sin(j*nt))    / (n^2*(j^2 - 1));
    sum2  = sum2 + (2*b(j)/j - cp(j)) * (cnt - cos(j*nt))      / (n^2*(j^2 - 1));
    sum3  = sum3 + (2*bp(j)/j + c(j)) * ((j*cnt - cos(j*nt)/j) / (n^2*(j^2 - 1)));
    sum4  = sum4 + (2*b(j)/j - cp(j)) * ((snt - sin(j*nt)/j)   / (n^2*(j^2 - 1)));
    sum9  = sum9 + (2*bp(j)/j + c(j)) * ((j*cnt - j*cos(j*nt))/(n*(j^2 - 1)));
    sum10 = sum10 + (2*b(j)/j - cp(j)) * ((snt - j*sin(j*nt))/(n*(j^2 - 1)));
    sum11 = sum11 + (2*bp(j)/j + c(j)) * (j*snt - sin(j*nt)) / (n*(j^2 - 1));
    sum12 = sum12 + (2*b(j)/j - cp(j)) * (cnt - cos(j*nt)) / (n*(j^2 - 1));
    sum23 = sum23 + u(j) * (j*snt - sin(j*nt)) / (n^2*(j^2-1));
    sum24 = sum24 + up(j) * (cnt - cos(j*nt)) / (n^2*(j^2-1));
    sum31 = sum31 + u(j) * (j*cnt - j*cos(j*nt)) / (n*(j^2-1));
    sum32 = sum32 + up(j) * (-snt + j*sin(j*nt)) / (n*(j^2-1));
end
for j = 2:6
    sum17 = sum17 + (2*mp(j)/j + r(j)) * (j*snt - sin(j*nt)) / (n^2*(j^2-1));
    sum18 = sum18 + (2*m(j)/j - rp(j)) * (cnt - cos(j*nt)) / (n^2*(j^2-1));
    sum19 = sum19 + (2*mp(j)/j + r(j)) * (-j*cnt + cos(j*nt)/j) / (n^2*(j^2-1));
    sum20 = sum20 + (2*m(j)/j - rp(j)) * (snt - sin(j*nt)/j) / (n^2*(j^2-1));
    sum25 = sum25 + (2*mp(j)/j + r(j)) * (j*cnt - j*cos(j*nt)) / (n*(j^2-1));
    sum26 = sum26 + (2*m(j)/j - rp(j)) * (-snt + j*sin(j*nt)) / (n*(j^2-1));
    sum27 = sum27 + (2*mp(j)/j + r(j)) * (j*snt - sin(j*nt)) / (n*(j^2-1));
    sum28 = sum28 + (2*m(j)/j - rp(j)) * (cnt - cos(j*nt)) / (n*(j^2-1));
end

xK_J2    = kp1*snt + kp2*cnt + (c0 + 2*n*K)*(1/n^2)*(1 - cnt) ...
           + 2*b0/n*(t - snt/n) + (2*bp(1) + c(1))*snt/(2*n^2) ...
           - (2.*bp(1) + c(1)).*t.*cnt./(2*n) - (2.*b(1) - cp(1)).*t.*snt./(2*n) ...
           + sum1 - sum2;
yK_J2    = 2*kp1*cnt - 2*kp2*snt - 2/n*(c0 + 2*n*K)*(t - snt/n) ...
           - 3*b0/2*t.^2 - 4*b0/n^2*cnt + 2*(2*bp(1) + c(1))*cnt/n^2 ...
           + (2*bp(1) + c(1))*t.*snt/n + (2*b(1) - cp(1))*(snt/n^2 - t.*cnt/n) ...
           + 2*sum3 + 2*sum4 - sum5 - sum6 + K*t + Kp;
zK_J2    = ep1*snt + ep2*cnt + d0/n^2*(1 - cnt) + d(1)/(2*n)*(snt/n - t.*cnt) + dp(1)/(2*n)*t.*snt ...
           + sum7 + sum8;
xdotK_J2 = n*kp1*cnt - n*kp2*snt + (c0 + 2*n*K)*snt/n + 2*b0/n*(1 - cnt) + (2*bp(1) + c(1))*t.*snt/2 ...
           - (2*b(1) - cp(1))*snt/(2*n) - (2*b(1) - cp(1))*t.*cnt/2 ...
           + sum9 + sum10;
ydotK_J2 = -2*n*kp1.*snt - 2*n*kp2.*cnt - 2/n*(c0 + 2*n*K)*(1 - cnt) ...
            - 4*b0*(t - snt/n) - (2*bp(1) + c(1))*snt/n ...
            + (2*bp(1) + c(1))*t.*cnt + (2*b(1) - cp(1))*t.*snt ...
            - 2*sum11 + 2*sum12 ...
            + b0*t - sum13 + sum14 + K;
zdotK_J2 = n*ep1.*cnt - n*ep2.*snt + d0/n.*snt + d(1)/n.*t.*snt + dp(1)/(2*n).*(snt + nt.*cnt) ...
           + sum15 + sum16;


xK_J3 = dx0./n.*snt + x0.*cnt + (r0 + 2.*n.*K1).*(1 - cnt)./n.^2 ...
        + 2.*m0./n.*(nt - snt)./n + (2.*mp(1) + r(1)).*snt./(2.*n.^2) ...
        - (2.*mp(1) + r(1)).*t.*cnt./(2.*n) - (2.*m(1) - rp(1)).*t.*snt./(2.*n) ...
        + sum17 - sum18;
yK_J3 = 2.*dx0./n.*cnt - 2.*x0.*snt - 2./n.*(r0 + 2.*n.*K1).*(t - snt./n) ...
        - 4.*m0.*(t.^2./2 + cnt./n.^2) - 2.*sum19 + 2.*sum20 ...
        + m0.*t.^2./2 - sum21 - sum22 + K1.*t + Kp1 ...
        + 2.*(2.*mp(1) + r(1)).*cnt./n.^2 + (2.*mp(1) + r(1)).*t.*snt./n ...
        + (2.*m(1) - rp(1)).*(snt./n.^2 - t.*cnt./n);
zK_J3 = dz0./n.*snt + z0.*cnt + u0./n.^2.*(1 - cnt) + u(1)./(2.*n).*(snt./n - t.*cnt) ...
        + up(1)./(2.*n).*t.*snt + sum23 + sum24;
xdotK_J3 = dx0.*cnt - n.*x0.*snt + (r0 + 2.*n.*K1).*snt./n + 2.*m0./n.*(1 - cnt) ...
           + (2.*mp(1) + r(1)).*t.*snt./2 - (2.*m(1) - rp(1)).*snt./(2.*n) ...
           - (2.*m(1) - rp(1)).*t.*cnt./2 + sum25 - sum26;
ydotK_J3 = -2.*dx0.*snt - 2.*n.*x0.*cnt - 2./n.*(r0 + 2.*n.*K1).*(1 - cnt) ...
           - 4.*m0.*(t - snt./n) - (2.*mp(1) + r(1)).*snt./n ...
           + (2.*mp(1) + r(1)).*t.*cnt + (2.*m(1) - rp(1)).*t.*snt ...
           - 2.*sum27 + 2.*sum28 + m0.*t - sum29 + sum30 + K1;
zdotK_J3 = dz0.*cnt - n.*z0.*snt + u0./n.*snt + u(1)./2.*t.*snt ...
           + up(1)./(2.*n).*(snt + nt.*cnt) + sum31 + sum32;

% Aplicar matriz de rotación
Rz_RAAN = [cos(Omega), -sin(Omega), 0;
           sin(Omega),  cos(Omega), 0;
                    0,           0, 1];
Rx_INC  = [1,        0,         0;
           0, cos(INC), -sin(INC);
           0, sin(INC),  cos(INC)];
xK    = zeros(length(t), 1);
yK    = zeros(length(t), 1);
zK    = zeros(length(t), 1);
xdotK = zeros(length(t), 1);
ydotK = zeros(length(t), 1);
zdotK = zeros(length(t), 1);

for j = 1:length(t)
    Rz_u = [cos(AOL(j)), -sin(AOL(j)), 0;
            sin(AOL(j)),  cos(AOL(j)), 0;
                      0,            0, 1];

    R_EUtoI = Rz_RAAN * Rx_INC * Rz_u;

    rJ2J3 = [   xK_J2(j) +    xK_J3(j);    yK_J2(j) +    yK_J3(j);    zK_J2(j) +    zK_J3(j)];
    vJ2J3 = [xdotK_J2(j) + xdotK_J3(j); ydotK_J2(j) + ydotK_J3(j); zdotK_J2(j) + zdotK_J3(j)];

    rJ2J3_I = R_EUtoI * rJ2J3;
    vJ2J3_I = R_EUtoI * vJ2J3;

    xK(j)    = rJ2J3_I(1);
    yK(j)    = rJ2J3_I(2);
    zK(j)    = rJ2J3_I(3);
    xdotK(j) = vJ2J3_I(1);
    ydotK(j) = vJ2J3_I(2);
    zdotK(j) = vJ2J3_I(3);
end

end