t = 0;

%% Define else where
phi=@(t) t*2; % roll %
theta=@(t) t; % pitch
si=@(t) t; % yaw
%%
%%Translational Conversion Matrix 
RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t)); 
      cos(si(t))*sin(phi(t))*sin(theta(t))-cos(si(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
      cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% body from global 

RGb= RbG' ;% global from body 

% Rotational Conversion Matrix 
%%
S=[ 1, 0, -sin(theta(t)); 
    0, cos(phi(t)), sin(phi(t))*cos(theta(t));
    0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
    % global from body 

S_inv= S^-1 ; % body from global 
%% Translational EOM 
m= 1; % mass in kg 

Fg=[0;0;m*9.81]; % airplane weight 
Fd=0; % drag =0 

% 4 Thrusts for each motor
T1=@(t) t*2 ;
T2=@(t) t*2 ;
T3=@(t) t*2 ;
T4=@(t) t*2 ;

F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 

Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM 


%% Rotational EOM 

I=[1 1 1; 
    1 1 1; 
    1 1 1]; %Moment of inertia 



Rot_G=0;
