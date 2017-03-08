% Shanmei Liu
% April 18, 2016

% if the initial value is 0, no periodic solution
% Seems the stable/ustable periodic solution deptend on the innitial value
%sol=dde23(@ddefun,21,1,[0 200]);
sol=dde23(@ddefun,21,0.9,[0 1200]);
plot(sol.x,sol.y);
xlabel('time t');
ylabel('N(t)');