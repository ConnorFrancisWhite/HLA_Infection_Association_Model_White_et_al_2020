%Function that will run the three pathogen ode model until a statiionary
%state with passed parameters.

function [t,y] = updated_model(i, fraction, M, deathRate, beta1, beta2, beta3, sigma1, sigma2, sigma3, mew1, mew2, mew3, m1, m2, m3);


%Default Parameters
if nargin == 0
    
    recovery = 0.001;
    beta = 8*(recovery+0.001);
    deathRate = 0.01;
    beta1 = beta;
    beta2 = beta;
    beta3 = beta;
    sigma1 = recovery;
    sigma2 = recovery;
    sigma3 = recovery;
    mew1 = recovery;
    mew2 = recovery;
    mew3 = recovery;
    m1 = 0;
    m2 = 0;
    m3 = 0;
    i = 0;
    M = 100;
end

%Initial Conditions

y = zeros(114,1);

%{
a = i/M;
b = fraction/M;
c = (1-a-b);
%}

a = i/M;
b = (1-a)*0.5;
c = (1-a)*0.5;


p1 = a;
p2 = b;
p3 = c;

P = [p1^2,2*p1*p2,p2^2,2*p2*p3,p3^2,2*p1*p3];

%Starting Infection
for k=1:1:6
   marker = (k-1)*19;
   y(1+marker) = 0.99*P(k);
   y(2+marker) = (1/300)*P(k);
   y(3+marker) = (1/300)*P(k);
   y(4+marker) = (1/300)*P(k);
    
end


time_steps = 200;
[t, y] = ode45(@ode_Three_Pathogens,[0 time_steps],[y],odeset('nonnegative',1:18,'AbsTol', 1e-14),[deathRate beta1 beta2 beta3 sigma1 sigma2 sigma3 mew1 mew2 mew3 m1 m2 m3]);

compareA = y(end,:);
compareB = y(1,:);
ll = 1;
while isequal(round(compareA,3,'significant'),round(compareB,3,'significant')) ~= 1
    ll = ll+1;
   
    ytemp = y(end,:);
    [x, ytemp] = ode45(@ode_Three_Pathogens,[t(end) t(end) + time_steps],[ytemp],odeset('nonnegative',1:18,'AbsTol', 1e-14),[deathRate beta1 beta2 beta3 sigma1 sigma2 sigma3 mew1 mew2 mew3 m1 m2 m3]);

    t = [t;x];
    y = [y;ytemp];
    
    
    compareA = y(end,:);
    compareB = y(end-100,:);
    for b = 1:1:length(y(end,:))
        
        if compareA(b) < 1*10^(-4)
           compareA(b) = 0;
        end
        if compareB(b) < 1*10^(-4)
           compareB(b) = 0;
        end
    end
    
   
    if ll == 100
        break;
    end
end




