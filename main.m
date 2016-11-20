clear all
clc
 
%%PARAMETERS
Ns = 100; 
d = zeros(8,Ns+1);
Ek = zeros(8, Ns+1);                                        % Initializing with all the values to 'zero'
u(1) = 17.3878;
u(2) = 0;
u(3) = 50;
Ek(:,1) = [12 0 80 0 0 0 0 0]';                           % Initial Condition
Uk(:,1) = [u(1) u(2) u(3) 0 0 0 0 0]';
d(:,1) = [0 0 0 0 0 0 0 0]';
Ts = 1;                                                  % Sampling Time 
Q = diag([0.01 0.00001 1 0.01 0.01 0.001 0.001 10]);     % Covariance matrix of model disturbances 
H = sqrt(Q);
I = eye(8);
C = [0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1];
R = diag([1 0.01 0.01 0.001 0.001 10]); 


%%TRUE STATES GENERATION
for j = 1:Ns+1
    randn('state', j)
    d(:,j) = sqrt(Q)*randn(8,1);
    
end

for i = 1:Ns 
    [T, e] = ode45(@(t,E) Model(t, E), [i*Ts,(i+1)*Ts], Ek(:, i)); % To solve the systems of equations
    Ek(:, i+1) = e(end, :)' + d(:,i+1);         % Adding gaussian noise to the states
end
Ek;

%%MEASUREMENT STATES
yk = zeros(6, Ns+1);

for k= 1:Ns+1
    vk(:,k) = sqrt(R)*randn(6, 1);
    yk(:,k) = C*Ek(:,k)+ vk(:,k);   
end


%Steady State Function is the Model Function%
fun = @steady_states_function;

%Linearisation%
A = Num_Jacobian(fun, Ek(:,1));
PHI = expm(A*Ts);

Gamma_D = eye(8);

Ek_hat = zeros(8, Ns+1);
Ek_hat_pred = zeros(8, Ns+1);
ek = zeros(6, Ns+1);
Pk_pred = zeros(8, 8, Ns+1);
Pk = zeros(8,8,Ns+1);

Ek_hat(:,1) = [12 0 80 0 0 0 0 0]'; 
ek(:,1) = yk(:,1) - C*Ek_hat(:,1);
Pk(:,:,1) = eye(8);


%% FOR THE PURSPOSE OF NEES ANALYSIS
alpha = 0.05;
Chi_1 = chi2inv(alpha, 8);
Chi_2 = chi2inv(1-alpha, 8);




%%FILTER IMPLEMENTATION

%% KALMAN FILTER
[Ek_kf ek_kf beta_kf P_kf Pp_kf est_err_kf RMSE_kf] = kf(Pk, ek)

%% EXTENDED KALMAN FILTER
[Ek_ekf ek_ekf beta_ekf P_ekf Pp_ekf est_err_ekf RMSE_ekf]= ekf(Pk, ek)

%% UNSCENTED KALMAN FILTER
[Ek_ukf ek_ukf beta_ukf P_ukf Pp_ukf est_err_ukf RMSE_ukf] = ukf(Pk)



%% RELEVANT GRAPHS
% figure

% plot(1:Ns+1, Ek(2,:), 'r'); hold on;
% plot(1:Ns+1, Ek_kf(2,:), 'b'); hold on;
% plot(1:Ns+1, Ek_ekf(2,:), 'm'); hold on;
% plot(1:Ns+1, Ek_ukf(2,:), 'c'); hold on;
% plot(1:Ns+1, yk(6,:), 'g')

%%Innovation Graphs
% figure

% plot(1:Ns+1, ek_kf(1,:), 'b'); hold on;
% plot(1:Ns+1, ek_ekf(1,:), 'm'); hold on;
% plot(1:Ns+1, ek_ukf(1,:), 'c'); hold on;


%% MEAN AND VARIANCE OF INNOVATION CALCULATIONS
for i =1:6
Mean_kf(i) = 0;
Mean_ekf(i) = 0;
Mean_ukf(i) = 0;

Var_kf(i) = 0;
Var_ekf(i) = 0;
Var_ukf(i) = 0;
end
for i = 1:6
    for j = 1:Ns+1
    
      Mean_kf(i) = Mean_kf(i) + ek_kf(i,j);
      Mean_ekf(i) = Mean_ekf(i) + ek_ekf(i,j);
      Mean_ukf(i) = Mean_ukf(i) + ek_ukf(i,j);
      
    end    
    Mean_kf(i) = Mean_kf(i)/(Ns+1)
    Mean_ekf(i) = Mean_ekf(i)/(Ns+1)
    Mean_ukf(i) = Mean_ukf(i)/(Ns+1)
    
    for k = 1:Ns+1
      
      Var_kf(i) = Var_kf(i) + (ek_kf(i,j)-Mean_kf(i))*(ek_kf(i,j)-Mean_kf(i));
      Var_ekf(i) = Var_ekf(i) + (ek_ekf(i,j)-Mean_ekf(i))*(ek_ekf(i,j)-Mean_ekf(i));
      Var_ukf(i) = Var_ukf(i) + (ek_ukf(i,j)-Mean_ukf(i))*(ek_ukf(i,j)-Mean_ukf(i));
        
    end 
    
    Var_kf(i) = Var_kf(i)/(Ns+1)
    Var_ekf(i) = Var_ekf(i)/(Ns+1)
    Var_ukf(i) = Var_ukf(i)/(Ns+1)
end
   
%%Spectral Radii Graphs

% plot(1:Ns+1, P_kf, 'b'); hold on;
% plot(1:Ns+1, P_ekf, 'm'); hold on;
% plot(1:Ns+1, P_ukf, 'g'); hold on;

plot(1:Ns+1, Pp_kf, 'b'); hold on;
plot(1:Ns+1, Pp_ekf, 'm'); hold on;
plot(1:Ns+1, Pp_ukf, 'g'); hold on;

%%Estimation Error
% figure

% plot(1:Ns+1, est_err_kf(8,:), 'b'); hold on;
% plot(1:Ns+1, est_err_ekf(8,:), 'm'); hold on;
% plot(1:Ns+1, est_err_ukf(8,:), 'c'); hold on;

% figure
% 
% plot(1:Ns+1, beta_kf, 'b'); hold on;
% plot(1:Ns+1, beta_ekf, 'm'); hold on;
% plot(1:Ns+1, beta_ukf, 'c'); hold on;
% plot(1:Ns+1, ones(Ns+1)*Chi_1, 'r'); hold on;
% plot(1:Ns+1, ones(Ns+1)*Chi_2, 'g'); hold on;

xlabel('Time')
ylabel('Spectral Radii of Updated Covariance')

legend('KF','EKF', 'UKF')

% RMSE_kf
% RMSE_ekf
% RMSE_ukf


