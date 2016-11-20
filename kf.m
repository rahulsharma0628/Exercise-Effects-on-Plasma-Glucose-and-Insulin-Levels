function [Ek_hat ek beta_kf Pu Pp est_err RMSE] = kf(Pk, ek )


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

yk = zeros(6, Ns+1);
% yk(:,1) = C*Ek(:,1)+ vk(:,1);

for i = 2:Ns+1 
    [T, e] = ode45(@(t,E) Model(t, E), [0,Ts], Ek(:, i-1)); % To solve the systems of equations
    Ek(:, i) = e(end, :)' + d(:,i);         % Adding gaussian noise to the states
    %d(:, i+1) = sqrt(Q)*randn(8,1);
%     yk(:,i) = C*Ek(:,i)+ vk(:,i);
end
Ek;

%%MEASUREMENT STATES
for k= 1:Ns+1
    vk(:,k) = sqrt(R)*randn(6, 1);
    yk(:,k) = C*Ek(:,k)+ vk(:,k);   
end


%%KALMAN FILTER IMPLEMENTATION
fun = @steady_states_function;
% E0 = [0 0 0 0 0 0 0 0];
% E_0 = fsolve(fun, E0)

Gamma_D = eye(8);

A = Num_Jacobian(fun, Ek(:,1));
PHI = expm(A*Ts);

Ek_hat = zeros(8, Ns+1);
Ek_hat_pred = zeros(8, Ns+1);
ek = zeros(6, Ns+1);
Pk_pred = zeros(8, 8, Ns+1);
Pk = zeros(8,8,Ns+1);

Ek_hat(:,1) = Ek(:,1);
ek(:,1) = yk(:,1) - C*Ek_hat(:,1);
Pk(:,:,1) = eye(8);
P = 0;
beta_kf(k) =(Ek(:,1) - Ek_hat(:,1))'*inv(Pk(:,:,1))*(Ek(:,1) - Ek_hat(:,1));
Pu(1) = trace(Pk(:,:,1));
Pp(1) = trace(Pk(:,:,1));
est_err(:,1) = Ek(:,1) - Ek_hat(:,1);

MSE = 0;
RMSE = 0;

for k = 2:Ns+1
    
    Ek_hat_pred(:,k) = PHI*Ek_hat(:,k-1);   
    Pk_pred(:,:,k) = PHI*Pk(:,:,k-1)*PHI' + Gamma_D*Q*Gamma_D';
    
    L = Pk_pred(:,:,k)*C'*(inv(C*Pk_pred(:,:,k)*C' + R));
    
    ek(:,k) = yk(:,k) - C*Ek_hat_pred(:,k);
    Ek_hat(:,k) = Ek_hat_pred(:,k) + L*ek(:,k);
    Pk(:,:,k) = (I - L*C)*Pk_pred(:,:,k);  
    
    beta_kf(k) =(Ek(:,k) - Ek_hat(:,k))'*inv(Pk(:,:,k))*(Ek(:,k) - Ek_hat(:,k));
    Pu(k) = max(abs(eig(Pk(:,:,k))));
    Pp(k) = max(abs(eig(Pk_pred(:,:,k))));
    
    est_err(:,k) = Ek(:,k) - Ek_hat(:,k);
    MSE = MSE + est_err(1,k)*est_err(1,k);
    
end
   RMSE = sqrt(MSE/Ns+1);

end
