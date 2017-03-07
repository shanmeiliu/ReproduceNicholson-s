% Shanmei Liu
% April 15, 2016

% Reproduce figures from paper below
% Hongying Shu, Lin Wang, Jianhong Wu, Global dynamics of Nicholson's blowflies equation revisited: Onset and termination of nonlinear oscillations, Journal of Differential Equations, Volume 255, Issue 9, 1 November 2013, Pages 2565-2586, ISSN 0022-0396, http://dx.doi.org/10.1016/j.jde.2013.06.020.
% (http://www.sciencedirect.com/science/article/pii/S0022039613002635


% parameters
gamma=1;
del=0.01;
p=15;
alpha=1;

tauspan = [0 75];
figure(1);
hold on;

for i=0:11
y0=pi/2+2*i*pi;
[t,y] = ode45( @ThetaN ,tauspan ,y0);
plot(t,y);

end
%disp([t,y])

b=gamma*(1-log(p/gamma)+del.*t);
s=t.*sqrt(b.^2-gamma^2);
plot(t,s,'r');

xlabel('\tau')
hold off