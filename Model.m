%%Defining the Model

function E_dot = Model(t, E)
I = E(1);
X = E(2);
G = E(3);
G_prod = E(4);
G_up = E(5);
I_e = E(6);
G_gly = E(7);
A = E(8);

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
u1 = 17.3878;
u2 = 0;
u3 = 50;
A_TH = -1.152*u3*u3 + 87.471*u3;

%Defining ODEs
E_dot(1)= -n*I + p4*u1 - I_e;
E_dot(2)= -p2*X + p3*(I - I_b);
E_dot(3)= -p1*(G - G_b) - X*G + W*(G_prod - G_gly - G_up)/Vol_G + u2/Vol_G;
E_dot(4)= a1*PVO2_max - a2*G_prod;
E_dot(5)= a3*PVO2_max - a4*G_up;
E_dot(6)= a5*PVO2_max - a6*I_e;
if (A >= A_TH) 
     E_dot(7) = k; 
else E_dot(7) = 0; 
end 
E_dot(8)= u3;

% To return a column vector
E_dot = E_dot';
end