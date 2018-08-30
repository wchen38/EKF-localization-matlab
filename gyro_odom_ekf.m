clear all
clc

imu_angular_vel = load('input_data/imu_angular_vel.csv');
imu_linear_accel = load('input_data/imu_accel.csv');
odom_velocity = load('input_data/odom_velocity.csv');
odom_pose = load('input_data/odom_pose.csv');
odom_pose_filtered = load('input_data/odom_pose_filtered.csv');
odom_theta = odom_pose(:,4);
odom_theta_filtered = odom_pose_filtered(:,4);

w = imu_angular_vel(:, 4);        %imu angular velocity in z axis
a = imu_linear_accel(:,2);        %imu linear acceleration in x axis

v = odom_velocity(:, 2);         %linear velocity from odom
odom_av = odom_velocity(:, 7);         %angular velocity from odom

START = 30;
END = 1015;
dt = 0.02;


%parameters
alpha1 = 0.1;
alpha2 = 0.2;
alpha3 = 0.5;
alpha4 = 0.6;

temp1 = [0.9^2, 0.9^2, 0.5^2, 0.01^2];
Qt = 0.05^2;%diag(temp1); 

X = [0;0;0; w(START)];
cov = [0.01^2 0.01^2 0.01^2 0.01^2];
P = diag(cov);
Pm = P;
Ht = [0 0 0 1];
      
Vt = [1 0;
      1 0;
      0 1;
      0 1]; 

 Xrec = [];

for k=START:END
    dt = (imu_angular_vel(k,1) - imu_angular_vel(k-1,1))/1000000000;
    theta = X(3);
    %theta = odom_theta_filtered(k);
    
    %-------------------------------prediction---------------------------
    Gt = [1 0 -v(k)*dt*sin(theta) 0;
          0 1 dt*cos(theta)       0;
          0 0 1                   dt;
          0 0 0                   1];
     
      Mt = [alpha1*v(k)^2 + alpha2*a(k)^2   0;
            0                               alpha3*v(k)^2 + alpha4*a(k)^2];
     
     
      Xm = X + [v(k)*cos(theta)*dt; v(k)*sin(theta)*dt; w(k)*dt; 0];
      Pm = Gt*P*Gt' + Vt*Mt*Vt';
      
      %-------------update-------------------------------
      innovation_cov = Ht*Pm*Ht' + Qt;
      Kt = Pm*Ht'/(innovation_cov);
      z = [w(k);];
      z_exp = [Xm(4)];
      innovation = z - z_exp;
      X = Xm + Kt*(innovation);
      X(3) = wrapToPi(X(3));
      P = (eye(4) - Kt*Ht)*Pm;
     % X = Xm;
      Xrec(:,k) = [X(1); X(2); X(3)];
end
plot(Xrec(1,:), Xrec(2,:)); hold on
plot(odom_pose_filtered(:,2), odom_pose_filtered(:,3));
plot( odom_pose(:,2), odom_pose(:,3) );
figure 

plot(Xrec(3,:)); hold on
plot(odom_pose_filtered(:,4));

