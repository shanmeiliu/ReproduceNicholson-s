% Shanmei Liu
% April 18, 2016

function yp = ddefun( t,y,Z )

ylag=Z;

tau=21;
delta=0.01;
p=15;
alpha=1;
gamma=1;

yp=exp(-delta*tau)*(p*ylag*exp(-alpha*ylag))-gamma*y;

end

