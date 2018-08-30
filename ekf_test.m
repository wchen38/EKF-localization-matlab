
clear all
clc

%------------parameters -----------------------------
deltaT = .02;
start = 30;
alphas = [.2 .03 .09 .08 0 0]; % motion model noise parameters

% measurement model noise parameters
sigma_range = .43;
sigma_bearing = .6;
sigma_id = 1;

Q_t = [sigma_range^2 0 0;
       0 sigma_bearing^2 0;
       0 0 sigma_id^2];

angle =[-1.57079637051:0.00436332309619:1.56643295288];
angle = angle';

landmarks = [0.372 1.85;
             3.95  1.95;
             5.33  5.73;
             2.23  -2.77;
             7     -1.95;
             8.36  1.08];
         
num_landmarks = size(landmarks, 1);

poseU = load('input_data/unfiltered_pose.txt');    %read in unfiltered pose data
poseU(:,1) = poseU(:,1)/1000000000; 
poseHat = load('input_data/filtered_pose.txt'); 
poseHat(:,1) = poseHat(:,1)/1000000000;
poseHat = poseHat(2:2514, :);

ranges = load('input_data/ranges.csv');            %read in the range of the scan
ranges(:,1) = ranges(:,1)/1000000000;                   %convert to seconds
ranges = ranges(2:2004, :)';

scan = ranges(2:721, 4);
linearVel = load('input_data/linearVel.txt');
angularVel = load('input_data/angularVel.txt');
quaternion1 = load('input_data/orientation.txt');

for i=1:size(quaternion1,1)
orientation(i,:) = quat2eul(quaternion1(i,2:5));
end


xc = scan.*cos(angle);
yc = scan.*sin(angle);
%plot(xc, yc);         %plot the scan data in cartesian coordiantes
%plot(landmarks(:,1), landmarks(:,2), 'o'); hold on
plot(poseHat(:,2),poseHat(:,3)); hold on
%inital guess of the robot pose
poseMean = [poseU(start, 2);
            poseU(start, 3);
            orientation(start,1)];
%inital guess of the pose covarience        
poseCov = zeros(3,3);
poseCov(:,:) = .01;

%--------------------algorithm ekf -----------------------
theta = orientation(start, 1);
for i = start:size(poseU, 1)
    if(i == 2004)
         break;
    end
    theta = poseMean(3); %I think this is wrong!!!
    deltaT = (poseU(i, 1) -  poseU(i-1, 1))/1000000000 ;
    
    u_t = [linearVel(i, 2); angularVel(i, 4)]; %linear and angular speed
    rot = deltaT * u_t(2);
    halfRot = rot / 2;
    trans = u_t(1) * deltaT;
    
    % calculate the movement Jacobian
    G_t = [1 0 trans * -sin(theta + halfRot);
           0 1 trans * cos(theta + halfRot);
           0 0 1];
    % calculate motion covariance in control space
    M_t = [(alphas(1) * abs(u_t(1)) + alphas(2) * abs(u_t(2)))^2 0;
           0 (alphas(3) * abs(u_t(1)) + alphas(4) * abs(u_t(2)))^2];
       
    % calculate Jacobian to transform motion covariance to state space
    V_t = [cos(theta + halfRot) -0.5 * sin(theta + halfRot);
           sin(theta + halfRot) 0.5 * cos(theta + halfRot);
           0 1];
       
    % calculate pose update from odometry
    poseUpdate = [trans * cos(theta + halfRot);
                  trans * sin(theta + halfRot);
                  rot];
              
    poseMeanBar = poseMean + poseUpdate;
    poseCovBar = G_t * poseCov * G_t' + V_t * M_t * V_t';
    
    A = ranges(1, :);               %range data time stamp
    timeStamp = poseU(i, 1);        %robot pose time stamp
    offset = [timeStamp - 0.04; timeStamp + 0.04];          %offset of pose time stamp
    
    currTimeStampIndex = find ( (A >= offset(1,1)) & (A <= offset(2,1)) );   %find the current range time stamp
    
    %currObsFeature = ranges(2:721, currTimeStampIndex(end));
    currObsFeature = ranges(2:721, i);
    currObsFeatureVec = [currObsFeature'; angle'; zeros(1,720)];
    %currObsFeatureVec = currObsFeatureVec';
    
    %xc = currObsFeature.*cos(angle);
    %yc = currObsFeature.*sin(angle);
  
    %plot(xc, yc, '-o');
    %pause(0.2);
    
    for j=start:size(currObsFeature)
           if(currObsFeature(j) > 2)
            continue;
           end
        % loop over all landmarks and compute MLE correspondence
            predZ = zeros(num_landmarks, 1, 3);
            predS = zeros(num_landmarks, 3, 3);
            predH  = zeros(num_landmarks, 3, 3);
            maxJ = 0;
            landmarkIndex = 1;
            
     
        for k=1:size(landmarks)
            xDist = landmarks(k, 1) - poseMeanBar(1);
            yDist = landmarks(k, 2) - poseMeanBar(2);
            q = xDist^2 + yDist^2;
            predZ(k,:,:) = [sqrt(q);
                           (atan2(yDist, xDist) - poseMeanBar(3));
                            0];
            predH(k,:,:) = [-xDist/sqrt(q) -yDist/sqrt(q) 0;
                                yDist/q        -xDist/q      -1;
                                0              0              0];
           predS(k,:,:) = squeeze(predH(k,:,:)) * poseCovBar ...
                               * squeeze(predH(k,:,:))' + Q_t;
           thisJ = det(2 * pi * squeeze(predS(k,:,:)))^(-0.5) * ...
                        exp(-0.5 * (currObsFeatureVec(:,j) - squeeze(predZ(k,:,:)))' ...
                        * inv(squeeze(predS(k,:,:))) ...
                        * (currObsFeatureVec(:,j) - squeeze(predZ(k,:,:))));
           if imag(thisJ) ~= 0
                    thisJ = 0;
           end
           
           if thisJ > maxJ
               maxJ = thisJ;
               landmarkIndex = k;
           elseif thisJ < 0
               disp("k less than 0");
               disp(thisJ);
           elseif thisJ == 0
               disp("k equals 0");
           end  
            
        end
        
        K = poseCovBar * squeeze(predH(landmarkIndex,:,:))' ...
                * inv(squeeze(predS(landmarkIndex,:,:)));
        % update pose mean and covariance estimates
        poseMeanBar = poseMeanBar + K * (currObsFeatureVec(:,j) - squeeze(predZ(landmarkIndex,:,:)));
        poseCovBar = (eye(3) - (K * squeeze(predH(landmarkIndex,:,:)))) * poseCovBar;
    end
    
    % update pose mean and covariance
    poseMean = poseMeanBar;
    poseMean(3) = poseMean(3);
    poseEstRec(:, i) = poseMean;
    poseCov = poseCovBar;
    poseCovRec(:, i) = [poseCov(1,1); poseCov(2,2); poseCov(3,3)];
    %theta = poseMean(3);
    
end





