%Function that will run the two pathogen ode model until a statiionary
%state with passed parameters.
function [t,y] = updated_model(i, M, deathRate, beta1, beta2, sigma1, sigma2, mew1, mew2, m1, m2, alpha, c);


%Default Parameters
if nargin == 0
    deathRate = 0.01;
    recovery = 0.02;
    beta = 0.06;
    sigma1 = recovery;
    sigma2 = recovery;
    
    mew1 = recovery;
    mew2 = recovery;
    beta1 = beta;
    beta2 = beta;
   
   
  
    m1 = 0;
    m2 = 0;
    
    %Allele Frequency for allele 1
    i = 60;
    M = 100;
   
    alpha = 0;
    c = 0.;

end


%Initial Conditions

y = zeros(27,1);

a = i/M;

b = (1-a);



p1 = a;
p2 = b;

P = [p1^2,2*p1*p2,p2^2];

%Starting Infection
for k=1:1:3
   marker = (k-1)*9;
   y(1+marker) = 0.99*P(k);
   y(2+marker) = (1/200)*P(k);
   y(3+marker) = (1/200)*P(k);
    
end


time_steps = 200;

[t, y] = ode45(@ode,[0 time_steps],[y],odeset('nonnegative',1:18,'AbsTol', 1e-14),[deathRate beta1 beta2 sigma1 sigma2 mew1 mew2 m1 m2, alpha, c]);


%This section makes sure the simulation reaches a stationary state.
compareA = y(end,:);
compareB = y(1,:);
ll = 1;
while isequal(round(compareA,3,'significant'),round(compareB,3,'significant')) ~= 1
    ll = ll+1;
    ytemp = y(end,:);
    [x,ytemp] = ode45(@ode,[t(end) t(end)+time_steps],[ytemp],odeset('nonnegative',1:18,'AbsTol', 1e-14),[deathRate beta1 beta2 sigma1 sigma2 mew1 mew2 m1 m2, alpha, c]);
    t = [t;x];
    y = [y;ytemp];
    
    compareA = y(end,:);
    compareB = y(end-100,:);
    for b = 1:1:length(y(end,:))
        
        if compareA(b) < 1*10^(-5)
           compareA(b) = 0;
        end
        if compareB(b) < 1*10^(-5)
           compareB(b) = 0;
        end
    end
    
    if ll == 100
        break;
    end
end



