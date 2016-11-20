%uncented kalman filter

function [Xukf innov_ukf beta_ukf P Pp est_err RMSE]  = ukf(Pk)

%%  PARAMETERS
Ns = 100; 
d = zeros(8,Ns+1);
Ek = zeros(8, Ns+1);                                        % Initializing w(i)th all the values to 'zero'
u(1) = 17.3878;
u(2) = 0;
u(3) = 50;
Ek(:,1) = [12 0 80 0 0 0 0 0]';                           % Initial Condition
d(:,1,1) = [0 0 0 0 0 0 0 0]';
Ts = 1;                                                  % Sampling Time 
Q = diag([0.01 0.00001 1 0.01 0.01 0.001 0.001 10]);     % Covariance matrix of model disturbances 
H = sqrt(Q);
I = eye(8);
C = [0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1];
R = diag([1 0.01 0.01 0.001 0.001 10]); 

%%TRUE STATES
for j = 1:Ns+1
    randn('state', j)
     d(:,j) = sqrt(Q)*randn(8,1);    
end

for i = 1:Ns 
    [T, e] = ode45(@(t,E) Model(t, E), [0,Ts], Ek(:, i)); % To solve the systems of equations
    Ek(:, i+1) = e(end, :)' + d(:,i+1);         % Adding gaussian noise to the states
    %d(:, i+1) = sqrt(Q)*randn(8,1);
end
%%MEASUREMENT STATES
yk = zeros(6, Ns+1);

for k= 1:Ns+1
    vk(:,k) = sqrt(R)*randn(6, 1);
    yk(:,k) = C*Ek(:,k)+ vk(:,k);   
end

%%UNSCENTED KALMAN FIlTER IMPLEMENTATION
NS = Ns+1;  
n = length(Ek(:,1));
m = length(yk(:,1));

Ek_hat = zeros(8, Ns+1);
Ek_hat_pred = zeros(8, Ns+1);
Pk_pred = zeros(8, 8, Ns+1);

Ek_hat(:,1) = Ek(:,1); 
ek(:,1) = yk(:,1) - C*Ek_hat(:,1);
Pk(:,:,1) = eye(8);

Xukf(:,1) = Ek(:,1);
%P = eye(n);
P_updated = Pk(:,:,1);
P(1) = max(abs(eig(P_updated)));
Pp(1) = max(abs(eig(P_updated)));
Ts = 1;

beta_ukf(1) = (Ek_hat(:,1) - Xukf(:,1))'/P_updated*(Ek_hat(:,1) - Xukf(:,1));

M = n + n + m;
K = 1;
rho = sqrt(M+K);

Pa = zeros(M,M);
Pa(n+1:n+n,n+1:n+n) = Q;
Pa(n+n+1:n+n+m,n+n+1:n+n+m) = R;

w = 1/2/(M+K)*ones(2*M+1,1);
w(1) = K/(M+K);
MSE = 0;
RMSE = 0;

for i = 2:NS
    
    Xa = [Xukf(:,i-1)' zeros(1,n) zeros(1,m)]';
    Pa(1:n,1:n) = P_updated;
    
    Xsam(1:M,1) = Xa;
    
    for j = 1:M
        Chi = zeros(M,1);
        Chi(j) = 1;
        
        Xsam(1:M,j+1) = Xa + rho*sqrtm(Pa)*Chi;
        Xsam(1:M,j+1+M) = Xa - rho*sqrtm(Pa)*Chi;
    end
    
    Xsam;
    
    for k = 1:2*M+1
     
     kT(i)=(i-1)*Ts;
     Initial_ode=Xsam(1:n,k);
     [T Y1] = ode45(@(t,x) Model(t,x),[(i-2)*Ts,(i-1)*Ts],Initial_ode);
     Xsam(1:n,k) = Y1(end,:)' + Xsam(n+1:2*n,k);
%    Xsam(1:n,k) = steady_states_function(Xsam(1:n,k));
    
    Ysam(1:m,k) = C*Xsam(1:n,k) + Xsam(2*n+1:end,k);
    
    end
    
    Xukf(:,i) = zeros(n,1);
    Yukf(:,i) = zeros(m,1);
    
    for k = 1:2*M+1
    Xukf(:,i) = Xukf(:,i) + w(k)*Xsam(1:n,k);
    Yukf(:,i) = Yukf(:,i) + w(k)*Ysam(:,k);
    end
    
    Pepseps = zeros(n,n);
    Pepse = zeros(n,m);
    Pee = zeros(m,m);
    
    for k = 1:2*M+1
        eps(1:n,k) = Xsam(1:n,k) - Xukf(:,i);
        e(1:m,k) = Ysam(1:m,k) - Yukf(:,i);
        
        Pepseps = Pepseps + w(k)*eps(1:n,k)*eps(1:n,k)';
        Pepse = Pepse + w(k)*eps(1:n,k)*e(1:m,k)';
        Pee = Pee + w(k)*e(1:m,k)*e(1:m,k)';
    end
    innov_ukf(:,1) = yk(:,1) - Yukf(:,1);
    est_err(:,1) = Ek(:,1) - Xukf(:,1);
    L = Pepse*inv(Pee);
    innov_ukf(:,i) = yk(:,i) - Yukf(:,i);
    Xukf(:,i) = Xukf(:,i) + L*(yk(:,i) - Yukf(:,i));
    P_updated = Pepseps - L*Pee*L';
    
   % P(i) = trace(P_updated);    
    P(i) = max(abs(eig(P_updated)));
    Pp(i) = max(abs(eig(Pepseps)));
    
    
    beta_ukf(i) = (Ek(:,i) - Xukf(:,i))'*inv(P_updated)*(Ek(:,i) - Xukf(:,i));
    
    est_err(:,i) = Ek(:,i) - Xukf(:,i);
    MSE = MSE + est_err(1,i)*est_err(1,i);
    
end
   RMSE = sqrt(MSE/Ns+1);
    
end
