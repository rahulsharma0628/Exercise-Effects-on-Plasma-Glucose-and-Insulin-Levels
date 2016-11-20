function [Ek_hat ek beta_ekf Pu Pp est_err RMSE] = ekf(Pk, ek)


%%PARAMETERS
Ns = 100; 
d = zeros(8,Ns+1);
Ek = zeros(8, Ns+1);                                        % Initializing with all the values to 'zero'
u(1) = 17.3878;
u(2) = 0;
u(3) = 50;
Ek(:,1) = [12 0 80 0 0 0 0 0]';                           % Initial Condition
d(:,1) = [0 0 0 0 0 0 0 0]';
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

%% EXTENDED KALMAN FILTER IMPLEMENTATION
%Gamma_D = integral(Fd, 0, Ts, 'ArrayValued', true)
Gamma_D = eye(8);

fun = @steady_states_function;

A(:,:,1) = Num_Jacobian(fun, Ek(:,1));
PHI(:,:,1) = expm(A(:,:,1)*Ts);

Ek_hat = zeros(8, Ns+1);
Ek_hat_pred = zeros(8, Ns+1);
Pk_pred = zeros(8, 8, Ns+1);

Ek_hat(:,1) = Ek(:,1); 
ek(:,1) = yk(:,1) - C*Ek_hat(:,1);
Pk(:,:,1) = eye(8);

beta_ekf(1) =(Ek(:,1) - Ek_hat(:,1))'*inv(Pk(:,:,1))*(Ek(:,1) - Ek_hat(:,1));
Pu(1) = max(abs(eig(Pk(:,:,1))));
Pp(1) = max(abs(eig(Pk(:,:,1))));
est_err(:,1) = Ek(:,1) - Ek_hat(:,1);

MSE = 0;
RMSE = 0;



for k = 2:Ns+1
    
    A(:,:,k) = Num_Jacobian(fun, Ek(:,k));
    PHI(:,:,k) = expm(A(:,:,k)*Ts);
    
    Ek_hat_pred(:,k) = PHI(:,:,k)*Ek_hat(:,k-1);
    Pk_pred(:,:,k) = PHI(:,:,k)*Pk(:,:,k-1)*PHI(:,:,k)' + Gamma_D*Q*Gamma_D';
    
    L = Pk_pred(:,:,k)*C'*(inv(C*Pk_pred(:,:,k)*C' + R));
    
    ek(:, k) = yk(:,k) - C*Ek_hat_pred(:,k);
    Ek_hat(:,k) = Ek_hat_pred(:,k) + L*ek(:,k);
    Pk(:,:,k) = (I - L*C)*Pk_pred(:,:,k);  
    
    [T, e] = ode45(@(t,E) Model(t, E), [0,Ts], Ek(:, k-1));
    Ek(:, k) = e(end, :)' + d(:,k);
    
    beta_ekf(k) =(Ek(:,k) - Ek_hat(:,k))'*inv(Pk(:,:,k))*(Ek(:,k) - Ek_hat(:,k));
    Pu(k) = max(abs(eig(Pk(:,:,k))));
    Pp(k) = max(abs(eig(Pk_pred(:,:,k))));
    
    est_err(:,k) = Ek(:,k) - Ek_hat(:,k);
    MSE = MSE + est_err(1,k)*est_err(1,k);
    
end
   RMSE = sqrt(MSE/Ns+1);


end
