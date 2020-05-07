%The function that describes the ODE for the two pathogen model.

function dydt = diff(t,y,pars)

d = pars(1);
beta1 = pars(2);
beta2 = pars(3);

sigma1 = pars(4);
sigma2 = pars(5);

mu1 = pars(6);
mu2 = pars(7);

m1 = pars(8);
m2 = pars(9);

alpha = pars(10);
c = pars(11);

N11 = sum(y(1:9));
N12 = sum(y(10:18));
N22 = sum(y(19:27));


NI1 = 0;
NI2 = 0;


for i = 1:1:3

	marker = (i-1)*9;
	
	NI1 = NI1 + y(2+marker) + y(5+marker) + y(8+marker);
	NI2 = NI2 + y(3+marker) + y(5+marker) + y(7+marker);

end

delta = [1 0 ; 1 1 ; 0 1];

p1 = N11 + N12/2; 
p2 = N22 + N12/2; 

total = p1 + p2;

p1 = p1/total;
p2 = p2/total;


p = [p1*p1,2*p1*p2,p2*p2];


lambda1 = beta1*NI1;
lambda2 = beta2*NI2;

lambda = lambda1 + lambda2;

dydt = [];

for i=1:1:3
    marker = (i-1)*9;
    dydt=[dydt; [      
    
    %1 SS
    p(i)*d - lambda*y(1+marker) + (1-delta(i,1))*sigma1*y(2+marker)*(1-alpha) + (1-delta(i,2))*sigma2*y(3+marker)*(1-alpha) - d*y(1+marker);
    %2 IS
    lambda1*y(1+marker) - lambda2*y(2+marker)*c - (1-delta(i,1))*sigma1*y(2+marker) + (1-delta(i,2))*sigma2*y(5+marker)*(1-alpha) - delta(i,1)*mu1*y(2+marker) - d*y(2+marker);
    %3 SI
    lambda2*y(1+marker) - lambda1*y(3+marker)*c - (1-delta(i,2))*sigma2*y(3+marker) + (1-delta(i,1))*sigma1*y(5+marker)*(1-alpha) - delta(i,2)*mu2*y(3+marker) - d*y(3+marker);
    %4 RS 
    delta(i,1)*mu1*y(2+marker)*(1-alpha) - lambda2*y(4+marker) + (1-delta(i,2))*sigma2*y(7+marker)*(1-alpha) - d*y(4+marker);
    %5 II
    lambda2*y(2+marker)*c + lambda1*y(3+marker)*c - (1-delta(i,1))*sigma1*y(5+marker) - (1-delta(i,2))*sigma2*y(5+marker) - delta(i,1)*mu1*y(5+marker) - delta(i,2)*mu2*y(5+marker) - d*y(5+marker);
    %6 SR
    delta(i,2)*mu2*y(3+marker)*(1-alpha) - lambda1*y(6+marker) + (1-delta(i,1))*sigma1*y(8+marker)*(1-alpha) - d*y(6+marker);
    %7 RI
    lambda2*y(4+marker) - (1-delta(i,2))*sigma2*y(7+marker) + delta(i,1)*mu1*y(5+marker)*(1-alpha) - delta(i,2)*mu2*y(7+marker) - d*y(7+marker);
    %8 IR
    lambda1*y(6+marker) - (1-delta(i,1))*sigma1*y(8+marker) + delta(i,2)*mu2*y(5+marker)*(1-alpha) - delta(i,1)*mu1*y(8+marker) - d*y(8+marker);
    %9 RR
    delta(i,2)*mu2*y(7+marker) + delta(i,1)*mu1*y(8+marker) + alpha*(delta(i,1)*mu1 + delta(i,2)*mu2 + (1-delta(i,1))*sigma1 + (1-delta(i,2))*sigma2)*y(5+marker) + alpha*((1-delta(i,1))*sigma1 + delta(i,1)*mu1)*y(2+marker) + alpha*((1-delta(i,2))*sigma2 + delta(i,2)*mu2)*y(3+marker) +  alpha*((1-delta(i,2))*sigma2)*y(7+marker) + alpha*((1-delta(i,1))*sigma1)*y(8+marker) - d*y(9+marker);
    ]];
    
end

