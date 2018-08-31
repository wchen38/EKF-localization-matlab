clear all
clc

imu_angular_vel = load('input_data/snake/imu_angular_vel_snake.csv');
imu_linear_accel = load('input_data/snake/imu_accel_snake.csv');
odom_velocity = load('input_data/snake/odom_velocity_snake.csv');
odom_pose = load('input_data/snake/odom_pose_snake.csv');
odom_pose_filtered = load('input_data/snake/odom_pose_filtered_snake.csv');
odom_theta = odom_pose(:,4);
odom_theta_filtered = odom_pose_filtered(:,4);

w = imu_angular_vel(:, 4);        %imu angular velocity in z axis
a = imu_linear_accel(:,2);        %imu linear acceleration in x axis

v = odom_velocity(:, 2);         %linear velocity from odom
odom_av = odom_velocity(:, 7);         %angular velocity from odom

START = 30;
%END = 1015;     %square data set
END = 1869;      %snake data set
dt = 0.02;


%parameters
alpha1 = 0.5;
alpha2 = 0.4;
alpha3 = 0.1;
alpha4 = 0.1;

%covariance Q, measurement noise (eq. 7.15)
temp1 = [0.9^2, 0.9^2, 0.5^2, 0.01^2];
Qt = diag(temp1); 

%initial guess of the pose
X = [0;0;0; v(START); w(START); a(START) ];

%initial guess of the covariance of the state vector
cov = [0.01^2 0.01^2 0.01^2 0.01^2 0.01^2 0.01^2];
P = diag(cov);
Pm = P;

%Ht is a Jacobian of the measurement model (eq. 7.14)
% map the a priori state x_{k | k-1} into the observed space which is 
%the measurement
Ht = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

%the motion noise to be mapped into state space (eq. 7.11)
Vt = [0 0;
      0 0;
      1 0;
      0 1;
      1 0;
      0 1]; 

 Xrec = []; %record all the pose inorder to plot later

for k=START:END
    dt = (imu_angular_vel(k,1) - imu_angular_vel(k-1,1))/1000000000;
    theta = X(3);
    
    %-------------------------------prediction---------------------------
    %The Jacobian from the motion model (eq. 7.8)
    Gt = [1 0 -v(k)*dt*sin(theta) dt*cos(theta) 0 0;
          0 1 dt*cos(theta)    dt*sin(theta) 0 0;
          0 0 1                0            dt 0;
          0 0 0                1             0 dt;
          0 0 0                0             1 0;
          0 0 0                0             0 1];
     
      %To derive the covariance of the additional motion noise nor(0, R)
      %we first determine the covariance matrix of the noise in control
      %space
      Mt = [alpha1*v(k)^2 + alpha2*a(k)^2   0;
            0                               alpha3*v(k)^2 + alpha4*a(k)^2];
     
      %motion model, unicycle model, to predict the pose
      Xm = X + [v(k)*cos(theta)*dt; v(k)*sin(theta)*dt; w(k)*dt; 0; 0; 0];
      %predict the covarence
      Pm = Gt*P*Gt' + Vt*Mt*Vt';
      
      %-------------update-------------------------------
      % innovation_cov: predict how much we should trust the measurement 
      % based on the a priori error covariance matrix P_{k | k-1} and 
      % the measurement covariance matrix R
      innovation_cov = Ht*Pm*Ht' + Qt; 
      
      %The Kalman gain is used to to indicate how much we trust the innovation
      Kt = Pm*Ht'/(innovation_cov);
      
      %measurement model
      z = [v(k); odom_av(k); w(k); a(k)];
      
      %expected measurements from our prediction
      z_exp = [Xm(4); Xm(5); Xm(5); Xm(6)];
      
      %innovation, difference between what we observe and what we expect
      innovation = z - z_exp;
      
      %update the pose
      X = Xm + Kt*(innovation);
      X(3) = wrapToPi(X(3));
      
      %update the covarence
      P = (eye(6) - Kt*Ht)*Pm;
      Xrec(:,k) = [X(1); X(2); X(3)];
end
plot(Xrec(1,:), Xrec(2,:)); hold on
plot(odom_pose_filtered(:,2), odom_pose_filtered(:,3));
plot( odom_pose(:,2), odom_pose(:,3) );
figure 

plot(Xrec(3,:)); hold on
plot(odom_pose_filtered(:,4));

