clc;
clear; 

% Orange prop is a 2 blade with 5in diameter and 4in Pitch
orangepropRPM = [2970 6066 9000 12000 15200 3040 6150 9000 12000 14000 3109 6072 9116 12063 15109];
orangepropThrust = [8 34 80 133.3 214.7 9 35.3 72.2 127.5 173.5 6.8 31.67 77.07 131.33 204];

% Black Prop is a 2 blade with 6in diameter and 3in pitch 
blackpropRPM = [3072 6098 8950 11974 12959 15000 2900 6010 9170 12000 14030 15000 3104 6005 9003 12008 14950];
blackpropThrust = [5.33 28 76 146 200 236 8.13 32.33 67.73 118.8 160 210 8.67 28.73 72 142 226]; 

% Green Prop is a three blade with 5in diameter and 4in pitch
greenpropRPM = [2097 5950 9000 12070 15400 16200 300 6130 9090 12030 14000 16000 3009 6145 8980 11915 15060 16191];
greenpropThrust = [6.7 27.3 67.3 118.7 193.3 214.66 9.13 32.3 67.7 118.8 160 210 6 25.4 64.27 114.67 181.6 207.73]; 

% used to find the power as a function of RPM to calculate torque
% T = 
rpm = [3009 6145 8980 11915 15060 16191 3109 6072 9116 12063 15109 3104 6005 9003 12008 14950];
electrical_power = [1.197 4.309 10.363 21.316 40.091 50.076 1.425 4.963 12.903 26.999 51.62 1.313 4.427 11.771 25.188 49.335];
radpersec = rpm * (pi/30);
torque = electrical_power./radpersec ;  

% Finding the coefficients for the second order equation for the props 
% it can be set up as in equation in the form below
% coeff_x(1,1)*w^2 + coeff_x(1,2)*w + coeff_x(1,1)
coeff_orange = polyfit (orangepropRPM, orangepropThrust, 2);
coeff_black = polyfit (blackpropRPM, blackpropThrust, 2);
coeff_green = polyfit (greenpropRPM, greenpropThrust,2);

coeff_torque = polyfit (radpersec, torque,1);

rpm = 0:0.01:25000; % set RPM range we want to extrapolate our data to 
radsec = rpm.*(pi/30); 

%corrected values are used to find kt
%assuming that coefficients for less than w^2 are 0 
orangeprop = (coeff_orange(1,1)*(rpm.^2) + (coeff_orange(1,2))*rpm - coeff_orange(1,3));
orangeprop_corrected = coeff_orange(1,1)*(rpm.^2); 

blackprop = (coeff_black(1,1)*(rpm.^2) + (coeff_black(1,2))*rpm - coeff_black(1,3));
blackprop_corrected = coeff_black(1,1)*(rpm.^2); 

greenprop = (coeff_green(1,1)*(rpm.^2) + (coeff_green(1,2))*rpm - coeff_green(1,3));
greenprop_corrected = coeff_green(1,1)*(rpm.^2); 

output_torque = coeff_torque(1,1)*radsec + coeff_torque(1,2);
%output_torque_corrected = coeff_torque(1,1)*radsec.^2;
%nour needs:
% Thrust=kt*w^2
% Torque=kd*w

kt = ["2 blade 5x4" coeff_orange(1,1) "g/rpm^2"; "2 blade 6x3" coeff_black(1,1) "g/rpm^2"; "3 blade 5x4" coeff_green(1,1) "g/rpm^2"] %kt values for each prop 


figure(1)
plot(orangepropRPM, orangepropThrust, "o")
hold on 
plot (rpm, orangeprop)
hold on 
plot(rpm, orangeprop_corrected)
hold on 
xlabel("Speed of Motor (Rev/min)")
ylabel("Thrust produced (g)")
title('2 Blade 5"x4" Prop')
grid on 


figure(2)
plot(blackpropRPM, blackpropThrust, "o")
hold on 
plot (rpm, blackprop)
hold on
plot (rpm,blackprop_corrected)
hold on 
xlabel("Speed of Motor (Rev/min)")
ylabel("Thrust produced (g)")
title('2 Blade 6"x3" Prop')
grid on 

figure(3)
plot(greenpropRPM, greenpropThrust, "o")
hold on 
plot (rpm, greenprop)
hold on
plot (rpm,greenprop_corrected)
hold on 
xlabel("Speed of Motor (Rev/min)")
ylabel("Thrust produced (g)")
title('3 Blade 5"x4" Prop')
grid on 

figure (4)
plot (radpersec,torque, "o")
hold on 
%plot (output_torque_corrected)
hold on 
plot (radsec, output_torque)
hold on 
xlabel("Angular Speed (rad/s) ")
ylabel("Output Torque (Nm)")
title('Torque as a function of Speed')
grid on




