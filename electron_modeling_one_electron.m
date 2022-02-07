global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²

boundary_x = 100e-9;
boundary_y = 200e-9;
n_particle = 10;
nSims = 2000;

mass_eff =  0.26*C.m_0;
T = 300;
Vth = sqrt((2* C.kb *T)/mass_eff);

tao_mn = 0.2e-12;
l_mfp = tao_mn*Vth;

t_step = 7e-15;

x_pos = (boundary_x)*rand;
y_pos = (boundary_y)*rand;
random_angle = (360)*rand;

Vx = Vth*cosd(random_angle);
Vy = Vth*sind(random_angle);

% if random_angle > 0 && random_angle <= (90)
%     Vx = Vth*cosd(random_angle);
%     Vy = Vth*sind(random_angle);
% elseif random_angle > (90) && random_angle <= 180
%     Vx = Vth*cosd(random_angle);
%     Vy = Vth*sind(random_angle);        
% elseif random_angle > 180 && random_angle <= (270)
%     Vx = Vth*cosd(random_angle);
%     Vy = Vth*sind(random_angle);        
% elseif random_angle > (270) && random_angle <= (360)
%     Vx = Vth*cosd(random_angle);
%     Vy = Vth*sind(random_angle);
% end

for n = 1: nSims
    
    prev_x_pos = x_pos;
    prev_y_pos = y_pos;
   
    Sx = Vx*t_step;
    Sy = Vy*t_step;
    
    x_pos = prev_x_pos + Sx;
    y_pos = prev_y_pos + Sy;
    
    if x_pos <= 0 
        x_pos = boundary_x;
        prev_x_pos = x_pos; 
    elseif x_pos >= boundary_x
        x_pos = 0;
        prev_x_pos = x_pos;
    end    
    if y_pos <= 0 
        y_pos = 0;
        Vy = Vy/(-1);
        
    elseif  y_pos >= boundary_y
        y_pos = boundary_y;
        Vy = Vy/(-1);
    end

    figure(1)
    plot([prev_x_pos x_pos],[prev_y_pos y_pos],'r');
    hold on
    xlim([0 100e-9]);
    ylim([0 200e-9]);
    pause(0.01);    

end
    
    

