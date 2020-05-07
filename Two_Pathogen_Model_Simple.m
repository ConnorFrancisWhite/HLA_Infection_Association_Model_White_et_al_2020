%Script to run the two pathogen model and produce a figure displaying how the 
%odds ratio changes with HLA frequency. Requires One_Run_Simple.m , ode.m, ORcalc.m 
%and Figure_Style.m  to run.


clearvars
close all;

M =100;


d = 0.01;
alpha = 0;
c = 0;

R0 =3;
recovery = 0.02;
beta_diff = 0.00;

sigma = recovery;
mu = recovery;
%beta = R0*(recovery+d);
beta = 0.06;
beta1 = beta + beta_diff;
beta2 = beta;

sigma1 = recovery;
sigma2 = recovery;

mew1 = mu;
mew2 = mu;

m1 = 0.0;
m2 = 0.0;

x = linspace(0,1,M+1);
x = x';

I1 = zeros(M+1,1);
I2 = zeros(M+1,1);
II = zeros(M+1,1);

R1 = zeros(M+1,1);
R2 = zeros(M+1,1);

OR1 = zeros(M+1,1);
OR2 = zeros(M+1,1);


for i=0:1:M
    
    [t,y] = One_Run_Simple(i,M,d, beta1, beta2, sigma1, sigma2, mew1, mew2, m1, m2, alpha,c);
    
    %Useful quantities
    NI1 = 0;
    NI2 = 0;
    NII = 0;
    NR1 = 0;
    NR2 = 0;

    for j = 1:1:3
        marker = (j-1)*9;
		
        NI1 = NI1 + y(:,2+marker) + y(:,5+marker) + y(:,8+marker);
        NI2 = NI2 + y(:,3+marker) + y(:,5+marker) + y(:,7+marker);
   
        
        NR1 = NR1 + y(:,4+marker) + y(:,7+marker) + y(:,9+marker);
        NR2 = NR2 + y(:,6+marker) + y(:,8+marker) + y(:,9+marker);
    end

    N11 = sum(y(:,1:9),2);
    N12 = sum(y(:,10:18),2);
    N22 = sum(y(:,19:27),2);
   
    N = [N11 N12 N22];
    total = N11 + N12 + N22;

    I1(i+1) = NI1(end); 
    I2(i+1) = NI2(end); 
    
    Z = [];
	temp = zeros(3,1);
	for g=1:1:3
		marker = (g-1)*9;
		Z = [Z;[y(length(y(:,1)),2+marker) + y(length(y(:,1)),3+marker) + y(length(y(:,1)),5+marker)  + y(length(y(:,1)),7+marker)  + y(length(y(:,1)),8+marker) ]];
	   
		temp(g) = N(1,g);
	end
	N = temp;
	
    
    %OR calculations
	PI1 = (Z(1)+Z(2))/(N(1) +  N(2));
    PI1N = (Z(3))/(N(3));
    PI2 = (Z(2)+Z(3))/(N(2) +  N(3));
    PI2N = (Z(1))/(N(1));    
    
    OR1(i+1) = ORcalc(PI1,PI1N);
    OR2(i+1) = ORcalc(PI2,PI2N);   

end        



one = zeros(M+1,1)+ 1;

figure( 'Position', [10 10 900 600]);
%Plot of equilibrium state
plot(x,I1,'k','LineWidth',3);
hold on 
plot(x,I2,'--k','LineWidth',3);

xlabel('p_{1}');
ylabel('Equilibrium Proportions');
ylim([0 1]);

L = {'I_{1}','I_{2}'};
legend(L);

Figure_Style;

%Plot of OR allele
figure( 'Position', [10 10 900 600]);


semilogy(x,OR1, 'k', 'LineWidth', 3);
hold on 

semilogy(x,OR2, '--k', 'LineWidth', 3);
hold on 
plot(x,one,'--k','LineWidth', 1);

xlabel('p_{1}');
ylabel('OR_{i}');

ylim([0.1,10]);
L = {'OR_{1}' 'OR_{2}'};
legend(L);
Figure_Style;










