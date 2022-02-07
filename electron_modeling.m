%% Part 1 Electron modeling
% Qiushi Chen 101049864

clear 
clc

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

boundary_x = 200e-9;
boundary_y = 100e-9;

n_particle = 1e4;
display_particle = 10;
nSims = 500;
T = 300;

prev_meanTemperature = T;


% Effective mass of electrons mn = 0.26*m0, m0 is the rest mass, m_0 = 9.10938215e-31  
mass_eff =  0.26*C.m_0;             % Effective mass of electrons mn
T = 300;                            % Temperature = 300K

% Question 1
% Thermal Velocity
Vth = sqrt((2* C.kb *T)/mass_eff); 
% Thermal Velocity = sqrt((2*Boltzmann constant * Temperature)/effective mass)
fprintf("Thermal Velocity = %d m/s \n", Vth);
% For a 2-D Model,thermal velocity has sqrt((2* C.kb *T)/mass_eff).
% For a 3-D Model,thermal velocity has sqrt((3* C.kb *T)/mass_eff).

x_pos = (boundary_x).*rand(1, display_particle);
y_pos = (boundary_y).*rand(1, display_particle);
random_angle = (360).*rand(1, display_particle);

% Question 2
% Mean Free Path
t_mn = 0.2e-12;                     % Mean time between collision 
meanFreePath = t_mn*Vth;            % Mean free path = Mean time between collision*velocity      
fprintf("Mean Free Path = %d m \n", meanFreePath);

% Pick a time_step of 5e-15 to have the special step to be smaller than
% 1/100 of region size
t_step = 7e-15;                     % time step
time = 0;                           % Initial Time = 0



Vx = Vth.*cosd(random_angle);
Vy = Vth.*sind(random_angle);
VTot = sqrt(Vx.^2 + Vy.^2);


figure(1)
subplot(211)
col = hsv(10);
hold on
xlim([0 200e-9]);
ylim([0 100e-9]);
xlabel('Length (m)');
ylabel('Width (m)');
for n = 1: nSims
    
    prev_time = time;
    time = prev_time + t_step;
    
    prev_x_pos = x_pos;
    prev_y_pos = y_pos;
    

    
    Sx = Vx.*t_step;
    Sy = Vy.*t_step;
    
    x_pos = prev_x_pos + Sx;
    y_pos = prev_y_pos + Sy;

    
    for k = 1:display_particle
        
        plot([prev_x_pos(k) x_pos(k)],[prev_y_pos(k) y_pos(k)],'color',col(k,:));
        
        % For the x direction use a periodic boundary condition where the particle jumps
        % to the opposite edge. i.e. if it reaches the right side it appears at the left with
        % the same velocity
        if x_pos(1,k) <= 0 
            x_pos(1,k) = boundary_x;
            y_pos(1,k) = prev_y_pos(1,k) ;
        elseif x_pos(1,k) >= boundary_x
            x_pos(1,k) = 0;
            y_pos(1,k) = prev_y_pos(1,k);
        end    
    
        
        % For the y direction use a boundary condition where the particle reflects at the
        % same angle (specular) and retains its velocity.
        if y_pos(1,k) <= 0 
            y_pos(1,k) = 0;
            Vy(1,k) = -(Vy(1,k));

        elseif  y_pos(1,k) >= boundary_y
            y_pos(1,k) = boundary_y;
            Vy(1,k) = -(Vy(1,k));
        end
        hold on
        

           
    end
    pause(0.001); 
    
    
    % Question 3 
    % Temperature plot
    % Calculate and display the semiconductor temperature on the plot at a fixed time
    % interval and verify that it stays constant.
    meanTemperature = mean(((VTot.^2).* mass_eff)./(2*C.kb));

    subplot(212)
    plot([prev_time time], [prev_meanTemperature meanTemperature],'r')
    ylim([100 500]);
    hold on
    
    prev_meanTemperature = meanTemperature;
    
    i=1;
    
    
end