%=====================================================
%      F16 nonlinear model trim cost function
%   for longitudinal motion, steady level flight
% (cost = sum of weighted squared state derivatives)
%
% Original authors: T. Keviczky, Richard S. Russell
% Date:   April 29, 2002, November 7th, 2002
%
%  File made suitable for trimming S-function dynamics
%
%  L. Sonneveldt May 2006
%
%
%=====================================================

function [cost, Xdot, xu, uu] = trim_fun(UX0)

global altitude velocity fi_flag_Simulink

% Implementing limits:
% thrust limits
if UX0(1) > 30*pi/180
    UX0(1) = 30*pi/180;
elseif UX0(1) < -30*pi/180
    UX0(1) = -30*pi/180;
end;

% elevator limits
if UX0(2) > 25*pi/180
    UX0(2) = 25*pi/180;
elseif UX0(2) < -25*pi/180
    UX0(2) = -25*pi/180;
end;

% angle of attack limits
if (fi_flag_Simulink == 0)
  if UX0(3) > 45*pi/180
    UX0(3) = 45*pi/180;
  elseif UX0(3) < -10*pi/180
    UX0(3) = -10*pi/180;
  end
elseif (fi_flag_Simulink == 1)
  if UX0(3) > 90*pi/180
    UX0(3) = 90*pi/180;
  elseif UX0(3) < -20*pi/180
    UX0(3) = -20*pi/180;
  end
end

%  Aileron limits
if UX0(4) > 21.5*pi/180
    UX0(4) = 21.5*pi/180;
elseif UX0(4) < -21.5*pi/180
    UX0(4) = -21.5*pi/180;
end;

% Rudder limits
if UX0(5) > 30*pi/180
    UX0(5) = 30*pi/180;
elseif UX0(5) < -30*pi/180
    UX0(5) = -30*pi/180;
end;

if (fi_flag_Simulink == 1)

    rho0 = 1.225;
    temp = 288.15-altitude*0.0065;
    rho = rho0*exp(-9.80665/287.05/temp*altitude);
    qbar = 0.5*rho*velocity^2;
    ps = rho/rho0*101325;
    
    dLEF = (1.38*UX0(3)*180/pi - 9.05*qbar/ps + 1.45)*pi/180;
    
elseif (fi_flag_Simulink == 0)
    dLEF = 0.0;
end

% dth limits
if UX0(6) > 1
    UX0(6) = 1;
elseif UX0(6) < 0
    UX0(6) = 0;
end;

% Verify that the calculated leading edge flap
% have not been violated.
if (dLEF > 25*pi/180)
    dLEF = 25*pi/180;
elseif (dLEF < 0)
    dLEF = 0;
end;

% Initialize the other variables
%
q0 = cos(UX0(3)/2);
q1 = 0;
q2 = sin(UX0(3)/2);
q3 = 0;
p = 0; 
q = 0; 
r = 0;
pow = tgear(UX0(6));  % taken from Ying Huo's model
   
tu = 0;
xu = [velocity UX0(1) UX0(3) q0 q1 q2 q3 p q r 0 0 -altitude pow]';
uu = [UX0(6) UX0(2) UX0(4) UX0(5) dLEF fi_flag_Simulink]';

dx = feval('F16_trim', tu, xu, uu, 'derivs');
 
Xdot = dx;

% Create weight function
weight = [  2            ...%Vt_dot
            10           ...%beta_dot
            10           ...%alpha_dot
            10           ...%q0_dot
            10           ...%q1_dot
            10           ...%q2_dot
            10           ...%q3_dot
            10           ...%p_dot
            10           ...%q_dot
            10           ...%r_dot
            0            ...%x_dot
            0            ...%y_dot
            5            ...%z_dot
            50            ...%pow_dot
            ];

cost = weight*(Xdot.*Xdot);