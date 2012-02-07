clc;
clear all;
close all;

%create the binary
mex F16_dyn.c

global altitude velocity fi_flag_Simulink alpha_trim V_trim h_trim pow_trim de_trim da_trim dr_trim 
global pitch_trim yaw_trim roll_trim init_x der_trim del_trim

fi_flag_Simulink = 1;
%create the binary

%% Trim aircraft to desired altitude and velocity
%%
altitude =6000;
Mach = 0.5;
c = sqrt(1.4*287.06*(288.15-0.0065*altitude));
velocity = Mach*c
simtime = 60*30;

%% Initial Conditions for trim routine.
%% The following values seem to trim to most
%% flight condition.  If the F16 does not trim
%% Change these values.
beta = 0;         % -
elevator = 0*pi/180;       % elevator, rad
alpha = 10*pi/180;           % AOA, rad
rudder = 0;         % rudder angle, rad
aileron = 0;         % aileron, rad
dth = 0.2;

% Initial Guess for free parameters
UX0 = [beta; elevator; alpha; aileron; rudder; dth];

% Initializing optimization options and running optimization:
OPTIONS = optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',5e+04,'MaxIter',1e+04);

iter = 1;
while iter == 1
    
    feval('F16_trim', [], [], [], 'lincompile');
    load_system('F16_trim');
   
    [UX,FVAL,EXITFLAG,OUTPUT] = fminsearch('trim_fun',UX0,OPTIONS);
   
    [cost, Xdot, xu, uu] = trim_fun(UX);
    
    disp('Trim Values and Cost:');
    disp(['cost   = ' num2str(cost)])
    disp(['dth    = ' num2str(uu(1)) ' -'])  
    disp(['elev   = ' num2str(uu(2)*180/pi) ' deg'])
    disp(['ail    = ' num2str(uu(3)*180/pi) ' deg'])
    disp(['rud    = ' num2str(uu(4)*180/pi) ' deg'])
    disp(['alpha  = ' num2str(xu(3)*180/pi) ' deg'])
    disp(['dLEF   = ' num2str(uu(5)*180/pi) ' deg'])
    disp(['Vel.   = ' num2str(xu(1)) ' m/s']) 
    disp(['pow    = ' num2str(xu(14)) ' %']) 
        iter = 0;
    feval('F16_trim', [], [], [], 'term');
    UX0 = UX;
end

alpha_trim = xu(3);
V_trim = xu(1);
h_trim = altitude;
pow_trim = xu(14);
de_trim = uu(2);
del_trim = uu(2)*180/pi;
der_trim = uu(2)*180/pi;
da_trim = uu(3)*180/pi;
dr_trim = uu(4)*180/pi;
dlef_trim = uu(5);
init_x = xu(1:14);
%init_x = [V_trim 0 alpha_trim, 0 alpha_trim 0, 0 0 0, 0 0 -h_trim, pow_trim];
pitch_trim = alpha_trim*180/pi;
yaw_trim =0;
roll_trim =0;
%%
%sim('F16error_fcs');

%%
%controller gains 
c11 = 0;c12 = 0.5;
c13 = 0.5;c21 =1;
c22 = 1;c23 = 1;

k11 = 0;k12 = 0.2;
k13 = 0.2;k21 = 0.5;
k22 = 0;k23 = 0;
%filter parameters 
mag_q = 50*pi/180;rate_q = 25*pi/180;zeta_q = 1.0; omega_q = 12.0;
mag_r = 15*pi/180;rate_r = 30*pi/180;zeta_r = 1.0; omega_r = 20.0;
mag_de = 25*pi/180;rate_de = 60*pi/180;zeta_de = 1.0; omega_de = 20.0;
mag_da = 21.5*pi/180;rate_da = 80*pi/180;zeta_da = 1.0; omega_da = 60.0;
mag_dr = 30*pi/180;rate_dr = 120*pi/180;zeta_dr = 1.0; omega_dr = 60.0;
%adaptive gains, dead zones and emod
lift_bound_u = 1*pi/180; lift_bound_l = 0000.01*pi/180; lift_gain_u = 0.1; lift_gain_l = 0.01; lift_emod = 0.001;
side_bound_u = 1*pi/180; side_bound_l = 000.01*pi/180; side_gain_u = 0.1; side_gain_l = 0.01; side_emod = 0.001;
l_bound_u = 5*pi/180; l_bound_l = 0000.01*pi/180; l_gain_u = 10.1; l_gain_l = .01; l_emod = 0.001;
m_bound_u = 5*pi/180; m_bound_l = 0000.01*pi/180; m_gain_u = 00100.1; m_gain_l = 005.1; m_emod = 0.001;
n_bound_u = 5*pi/180; n_bound_l = 00000.01*pi/180; n_gain_u = 0.004; n_gain_l = 0.001; n_emod = 0.001;
%initital theta
%theta_ini = thetahat(end,:);
theta_ini = 0.0*ones(1,105000);
time_failure = 200;
size_failure = 0.125*21.5;  %0, 0.125, 0.25, 0.5
sign_failure = -1;

sim('F16_fcs');