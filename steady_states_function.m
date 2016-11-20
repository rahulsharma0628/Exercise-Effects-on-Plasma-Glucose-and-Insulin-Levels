function F = steady_states_function(x)


%Values of all the Parameters 
p1 = 0.035;
p2 = 0.05;
p3 = 0.000028;
p4 = 0.098;
n = 0.142;
Vol_G = 117;
G_b = 80;
I_b = 12;
a1 = 0.00158;
a2 = 0.056;
a3 = 0.00195;
a4 = 0.0485;
a5 = 0.00125;
a6 = 0.075;
k = 0.0108;
T1 = 6;
PVO2_max = 40;
W = 60;
N = 60;
Ts = 1;

%Input Values
u(1) = 17.3878;
u(2) = 0;
u(3) = 50;
A_TH = -1.152*u(3)*u(3) + 87.471*u(3);

F(1)= -n*x(1)+ p4*u(1) - x(6);
F(2)= -p2*x(2) + p3*(x(1)- I_b);
F(3)= -p1*(x(3)- G_b) - x(2)*x(3)+ W*(x(4) - 0 - x(5))/Vol_G + u(2)/Vol_G ;
F(4)= a1*PVO2_max - a2*x(4);
F(5)= a3*PVO2_max - a4*x(5);
F(6)= a5*PVO2_max - a6*x(6);
if (x(8) >= A_TH) 
     F(7) = k; 
else F(7) = 0; 
end 
F(8)= u(3);

end