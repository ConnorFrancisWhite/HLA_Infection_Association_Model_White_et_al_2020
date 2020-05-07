%The function that describes the ODE for the three pathogen model.
function dydt = diff(t,y,pars)

d = pars(1);
beta1 = pars(2);
beta2 = pars(3);
beta3 = pars(4);
sigma1 = pars(5);
sigma2 = pars(6);
sigma3 = pars(7);
mu1 = pars(8);
mu2 = pars(9);
mu3 = pars(10);
m1 = pars(11);
m2 = pars(12);
m3 = pars(13);



N11 = sum(y(1:19));
N12 = sum(y(20:38));
N22 = sum(y(39:57));
N23 = sum(y(58:76));
N33 = sum(y(77:95));
N13 = sum(y(96:114));

NI1 = 0;
NI2 = 0;
NI3 = 0;

for i = 1:1:6
    marker = (i-1)*19;
    
    NI1 = NI1 + y(2+marker) + y(10+marker) + y(12+marker) + y(19+marker);
    NI2 = NI2 + y(3+marker) + y(8+marker) + y(13+marker) + y(18+marker);
    NI3 = NI3 + y(4+marker) + y(9+marker) + y(11+marker) + y(17+marker);
    
end


delta = [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1];




p1 = N11 + N12/2 + N13/2;
p2 = N22 + N12/2 + N23/2;
p3 = N33 + N13/2 + N23/2;

total = p1 + p2 + p3;

p1 = p1/total;
p2 = p2/total;
p3 = p3/total;


p = [p1*p1,2*p1*p2,p2*p2,2*p2*p3,p3*p3,2*p1*p3];


lambda1 = beta1*NI1;
lambda2 = beta2*NI2;
lambda3 = beta3*NI3;

lambda = lambda1 + lambda2 + lambda3;


dydt = [];


for i=1:1:6
    marker = (i-1)*19;
    dydt=[dydt; [      
    
    %1 SSS
    p(i)*(d + NI1*m1 + NI2*m2 + NI3*m3 )  - lambda*y(1+marker) + (1-delta(i,1))*sigma1*y(2+marker) + (1-delta(i,2))*sigma2*y(3+marker) + (1-delta(i,3))*sigma3*y(4+marker) - d*y(1+marker);
    %2 ISS
    lambda1*y(1+marker) - (1-delta(i,1))*sigma1*y(2+marker) - delta(i,1)*mu1*y(2+marker) - d*y(2+marker) - (1-delta(i,1))*y(2+marker)*m1;
    %3 SIS
    lambda2*y(1+marker) - (1-delta(i,2))*sigma2*y(3+marker) - delta(i,2)*mu2*y(3+marker) - d*y(3+marker) - (1-delta(i,2))*y(3+marker)*m2;
    %4 SSI
    lambda3*y(1+marker) - (1-delta(i,3))*sigma3*y(4+marker) - delta(i,3)*mu3*y(4+marker) - d*y(4+marker) - (1-delta(i,3))*y(4+marker)*m3;
    %5 RSS
    delta(i,1)*mu1*y(2+marker) - lambda2*y(5+marker) - lambda3*y(5+marker) + (1-delta(i,2))*sigma2*y(8+marker) + (1-delta(i,3))*sigma3*y(9+marker) - d*y(5+marker);
    %6 SRS
    delta(i,2)*mu2*y(3+marker) - lambda1*y(6+marker) - lambda3*y(6+marker) + (1-delta(i,1))*sigma1*y(10+marker) + (1-delta(i,3))*sigma3*y(11+marker) - d*y(6+marker);
    %7 SSR
    delta(i,3)*mu3*y(4+marker) - lambda1*y(7+marker) - lambda2*y(7+marker) + (1-delta(i,1))*sigma1*y(12+marker) + (1-delta(i,2))*sigma2*y(13+marker) - d*y(7+marker);
    %8 RIS
    lambda2*y(5+marker) - (1-delta(i,2))*sigma2*y(8+marker) - delta(i,2)*mu2*y(8+marker) - d*y(8+marker) - (1-delta(i,2))*y(8+marker)*m2;
    %9 RSI
    lambda3*y(5+marker) - (1-delta(i,3))*sigma3*y(9+marker) - delta(i,3)*mu3*y(9+marker) - d*y(9+marker) - (1-delta(i,3))*y(9+marker)*m3;
    %10 IRS
    lambda1*y(6+marker) - (1-delta(i,1))*sigma1*y(10+marker) - delta(i,1)*mu1*y(10+marker) - d*y(10+marker) - (1-delta(i,1))*y(10+marker)*m1;
    %11 SRI
    lambda3*y(6+marker) - (1-delta(i,3))*sigma3*y(11+marker) - delta(i,3)*mu3*y(11+marker) - d*y(11+marker) - (1-delta(i,3))*y(11+marker)*m3;
    %12 ISR
    lambda1*y(7+marker) - (1-delta(i,1))*sigma1*y(12+marker) - delta(i,1)*mu1*y(12+marker) - d*y(12+marker) - (1-delta(i,1))*y(12+marker)*m1;
    %13 SIR
    lambda2*y(7+marker) - (1-delta(i,2))*sigma2*y(13+marker) - delta(i,2)*mu2*y(13+marker) - d*y(13+marker) - (1-delta(i,2))*y(13+marker)*m2;
    %14 RRS
    delta(i,2)*mu2*y(8+marker) + delta(i,1)*mu1*y(10+marker) - lambda3*y(14+marker) + (1-delta(i,3))*sigma3*y(17+marker) - d*y(14+marker);
    %15 RSR
    delta(i,3)*mu3*y(9+marker) + delta(i,1)*mu1*y(12+marker) - lambda2*y(15+marker) + (1-delta(i,2))*sigma2*y(18+marker) - d*y(15+marker);
    %16 RRS
    delta(i,3)*mu3*y(11+marker) + delta(i,2)*mu2*y(13+marker) - lambda1*y(16+marker) + (1-delta(i,1))*sigma1*y(19+marker) - d*y(16+marker);
    %17 RRI
    lambda3*y(14+marker) - (1-delta(i,3))*sigma3*y(17+marker) - d*y(17+marker) - (1-delta(i,3))*y(17+marker)*m3;
    %18 RIR
    lambda2*y(15+marker) - (1-delta(i,2))*sigma2*y(18+marker) - d*y(18+marker) - (1-delta(i,2))*y(18+marker)*m2;
    %19 IRR
    lambda1*y(16+marker) - (1-delta(i,1))*sigma1*y(19+marker) - d*y(19+marker) - (1-delta(i,1))*y(19+marker)*m1;
    ]];

end

