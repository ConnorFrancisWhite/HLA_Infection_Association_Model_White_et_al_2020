%Script to run the three pathogen model and produce a figure displaying how the 
%odds ratio changes with HLA frequency. Requires One_Run_Three_Pathogen_Simple.m 
%, ode_Three_Pathogens.m, ORcalc.m and Figure_Style.m  to run.


clearvars
close all
deathRate = 0.01;


R0 = 4;
recovery = 0.02;
mu = 0.02;
%beta = R0*(recovery+deathRate);
beta = 0.06;
beta1 = beta;
beta2 = beta;
beta3 = beta;

sigma1 = recovery;
sigma2 = recovery;
sigma3 = recovery;
mew1 = mu;
mew2 = mu;
mew3 = mu;
m1 = 0.0;
m2 = 0.0;
m3 = 0.0;

M = 100;
x = linspace(0,1,M+1);

R1 = zeros(M+1,1);
R2 = zeros(M+1,1);
R3 = zeros(M+1,1);

I1 = zeros(M+1,1);
I2 = zeros(M+1,1);
I3 = zeros(M+1,1);

OR1 = zeros(M+1,1);
OR2 = zeros(M+1,1);
ORE = zeros(M+1,1);

fraction = 2;



for i=0:1:M
    p1 = i/M;
    p2 = (1 - p1)/2;
    p3 = (1 - p1)/2;
 
  
    [t,y] = One_Run_Three_Pathogens_Simple(i,fraction,M,deathRate, beta1, beta2, beta3, sigma1, sigma2,  sigma2, mew1, mew2, mew3, m1, m2, m3);

    %Useful quantities
    NI1 = 0;
    NI2 = 0;
    NI3 = 0;
    NR1 = 0;
    NR2 = 0;
    NR3 = 0;

    for j = 1:1:6
        marker = (j-1)*19;

        NI1 = NI1 + y(:,2+marker) + y(:,10+marker) + y(:,12+marker) + y(:,19+marker);
        NI2 = NI2 + y(:,3+marker) + y(:,8+marker) + y(:,13+marker) + y(:,18+marker);
        NI3 = NI3 + y(:,4+marker) + y(:,9+marker) + y(:,11+marker) + y(:,17+marker);
        NR1 = NR1 + y(:,5+marker) + y(:,8+marker) + y(:,9+marker) + y(:,14+marker) + y(:,15+marker) + y(:,17+marker) + y(:,18+marker);
        NR2 = NR2 + y(:,6+marker) + y(:,10+marker) + y(:,11+marker) + y(:,14+marker) + y(:,16+marker) + y(:,17+marker) + y(:,19+marker);
        NR3 = NR3 + y(:,7+marker) + y(:,12+marker) + y(:,13+marker) + y(:,15+marker) + y(:,16+marker) + y(:,18+marker) + y(:,19+marker);
    end

    N11 = sum(y(:,1:19),2);
    N12 = sum(y(:,20:38),2);
    N22 = sum(y(:,39:57),2);
    N23 = sum(y(:,58:76),2);
    N33 = sum(y(:,77:95),2);
    N13 = sum(y(:,96:114),2);

    N = [N11 N12 N22 N23 N33 N13];
    total = N11 + N12 + N22 + N23 + N33 + N13;

    I1(i+1) = NI1(end); 
    I2(i+1) = NI2(end); 
    I3(i+1) = NI3(end); 
    

    Z = [];
    temp = zeros(6,1);
    for g=1:1:6
        marker = (g-1)*19;
        Z = [Z;[y(length(y(:,1)),2+marker) + y(length(y(:,1)),3+marker) + y(length(y(:,1)),4+marker)  + y(length(y(:,1)),8+marker)  + y(length(y(:,1)),9+marker) + y(length(y(:,1)),10+marker) + y(length(y(:,1)),11+marker) + y(length(y(:,1)),12+marker) + y(length(y(:,1)),13+marker)  + y(length(y(:,1)),17+marker) + y(length(y(:,1)),18+marker) + y(length(y(:,1)),19+marker)]];

        temp(g) = N(1,g);
    end
    N = temp;
	
	
	%OR calculation
    
    PI1 = (Z(1)+Z(2)+Z(6))/(sum(Z));
    PI1N = (N(1)+N(2)+N(6)-Z(1)-Z(2)-Z(6))/(sum(N) - sum(Z));
    
    PI2 = (Z(2)+Z(3)+Z(4))/(sum(Z));
    PI2N = (N(2)+N(3)+N(4)-Z(2)-Z(3)-Z(4))/(sum(N) - sum(Z));
    
    PIE = sum(Z(2:6))/(sum(Z));
    PIEN = (sum(N(2:6)) - sum(Z(2:6)))/(sum(N) - sum(Z));
    
	OR1(i+1) = ORcalc(PI1,PI1N);
    OR2(i+1) = ORcalc(PI2,PI2N);
	ORE(i+1) = ORcalc(PIE,PIEN);
	
	
   
    
end
x = linspace(0,1,M+1);
x = x';
one = zeros(M+1,1)+ 1;


figure
semilogy(x,OR1, 'LineWidth', 3);
hold on;
semilogy(x,OR2, 'LineWidth', 3);
hold on;
semilogy(x,ORE, 'LineWidth', 3);
hold on;
plot(x,one,'--k','LineWidth', 0.5);
xlabel('p_{1}');
ylabel('OR_{i}');
ylim([0.1 10]);
L = {'OR_{1}','OR_{2}','OR_{2,3}'};
legend(L);
Figure_Style;

figure( 'Position', [10 10 900 600]);
semilogy(x,OR1,'k', 'LineWidth', 3);
hold on;
plot(x,one,'--k','LineWidth', 0.5);
xlabel('p_{1}');
ylabel('OR_{i}');
ylim([0.1 10]);
L = {'OR_{1}'};
legend(L);
Figure_Style;


figure;
%Plot of equilibrium state
plot(x,I1,'LineWidth',3);
hold on 
plot(x,I2,'LineWidth',3);
hold on
plot(x,I3,'LineWidth',3);




xlabel('p_{1}');
ylabel('I_{i}');
L = {'I_{1}','I_{2}','I_{3}'};
legend(L);
Figure_Style;



