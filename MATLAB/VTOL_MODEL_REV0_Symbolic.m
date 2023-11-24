close all 
clear all 


%% Define else where
syms phi(t) theta(t) si(t) t 
phi(t); % phi(t) roll %
theta(t); % pitch
si(t); % yaw



%%
%%Translational Conversion Matrix 
RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t)); 
      cos(si(t))*sin(phi(t))*sin(theta(t))-cos(si(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
      cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% body from global 

RGb= RbG' ;% global from body 

% Rotational Conversion Matrix 
%%
% Rotational Conversion Matrix 
S=[ 1, 0, -sin(theta(t)); 
    0, cos(phi(t)), sin(phi(t))*cos(theta(t));
    0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
    % global from body 

S_inv= S^-1 ; % body from global 
%% Translational EOM 
syms m g T1(t) T2(t) T3(t) T4(t) W1(t) W2(t) W3(t) W4(t) 
m; % mass in kg  

Fg=[0;0;m*g]; % airplane weight 
Fd=0; % drag =0 

% 4 speeds for each motor
W1(t) ;
W2(t)  ;
W3(t)  ;
W4(t)  ;

% Thrust as a function of time given motor speed
T1(t) ;
T2(t) ;
T3(t) ;
T4(t) ;

F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 

Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM 


%% Rotational EOM 

syms  Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz 
syms P1x P1y P1z P2x P2y P2z P3x P3y P3z P4x P4y P4z

I=[Ixx Ixy Ixz; 
    Ixy Iyy Iyz; 
    Ixz Iyz Izz]; %Moment of inertia 

% Position of each motor relative to body frame from CG
P1=[P1x;P1y;P1z];
P2=[P2x;P2y;P2z];
P3= [P3x;P3y;P3z];
P4=[P4x;P4y;P4z];

%Motors moment arm contribution in body axis frame
Tm1=cross(P1,[0;0;-T1(t)]);
Tm2=cross(P2,[0;0;-T2(t)]);
Tm3=cross(P3,[0;0;-T3(t)]);
Tm4=cross(P4,[0;0;-T4(t)]);

Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;

%kd is drag torque proportionality constant
kd=1;
% rotation direction of each motor. +ve=cw, -ve=ccw

syms M1 M2 M3 M4


Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];

Tm_total= Tm_pitch_roll + Tm_yaw ;

%gyroscopic effect: 


 % torque moment of inertia 
syms phi_derivative(t) theta_derivative(t) si_derivative(t) J_r

T_g=[J_r*theta_derivative(t)*(M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t); 
    -1*J_r*phi_derivative(t)*(M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t);
    0];
T_g=0;
w=[phi_derivative(t);theta_derivative(t);si_derivative(t)];

% third term:
T_g_2=cross(w,I*w);

w_dot_symb = I\((Tm_total)-(T_g)-(T_g_2));
w_dot = convert_symbols_to_floats(w_dot_symb);


theta_dot_der = w_dot(1);
theta_der = theta_derivative(float_t);
phi_dot_der = w_dot(2);
phi_der = phi_derivative(float_t);
si_dot_der = w_dot(3);
si_der = si_derivative(float_t);

solv_array = [theta_dot_der;theta_der;phi_dor_der;phi_der;si_dot_der;si_der];
[theta_dot,theta,phi_dot,phi,si_dot,si] = ode45(solv_array,[0 10],[0 0 0 0 0 0])

%% Function to convert symbols to floats 
function result_equation = convert_symbols_to_floats(input_equation)
    % Extract all the symbols from the input equation
    symbol_list = symvar(input_equation);
    
    % Create a cell array to store original symbols and undefined floats
    symbol_cells = cell(1, numel(symbol_list)*2);
    
    % Replace each symbol with an undefined float and store the mapping
    for i = 1:numel(symbol_list)
        symbol_cells{i} = symbol_list(i);
        symbol_cells{i+numel(symbol_list)} = sym(['float_', char(symbol_list(i))]);
        input_equation = subs(input_equation, symbol_cells{i}, symbol_cells{i+numel(symbol_list)});
    end
    
    % Return the equation with undefined floats
    result_equation = input_equation;
end