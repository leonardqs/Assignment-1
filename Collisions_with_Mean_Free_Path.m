%% Part 2 Collisions with Mean Free Path
% Qiushi Chen 101049864


global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      %metres (32.1740 ft) per s²

boundary_x = 200e-9;                % Boundary length
boundary_y = 100e-9;                % Boundary width
t_step = 7e-15;                     % Time Step
n_particle = 1e4;
display_particle = 10;              % Display Particle in the box
nSims = 1000;                       % Simulation Time
T = 300;                            % Default Temperature
t_mn = 0.2e-12;                     % Mean time between collision 
prev_meanTemperature = T;


mass_eff =  0.26*C.m_0;
initial_Vth = sqrt((2* C.kb *T)/mass_eff);

% Generate the initial electron position
x_pos = (boundary_x).*rand(1, display_particle);    % Initial x position
y_pos = (boundary_y).*rand(1, display_particle);    % Initlal Y position     


% Question 1 
% Plot the Maxwell-Boltzmann distribution for each electron velocity component in a histogram

% Generate the initial electron thermal velocity histogram
Vx_hist = initial_Vth.*randn(n_particle,1);       % Initial x velocity
Vy_hist = initial_Vth.*randn(n_particle,1);       % Initial y velocity
Vth_hist = sqrt(Vx.^2 + Vy.^2);                   % Thermal velocity
figure(1)
hist(Vth_hist,100);                               % Histogram
xlabel('Thermal Velocity (m/s)');
ylabel('Relative number of molecule per velocity');


Vx = (initial_Vth/sqrt(2)).*randn(1,display_particle);
Vy = (initial_Vth/sqrt(2)).*randn(1,display_particle);
VTot = sqrt(Vx.^2 + Vy.^2);
time = 0;

% exponential scattering probability
Pscat = 1- exp(-(t_step/t_mn));

figure(2)
col = hsv(10);
hold on
xlim([0 200e-9]);
ylim([0 100e-9]);

xlabel('Length (m)');
ylabel('Width (m)');
for n = 1: nSims
    
    % Update the electron position and velocity each time step
    prev_time = time;
    time = prev_time + t_step;
    
    prev_x_pos = x_pos;
    prev_y_pos = y_pos;

    
    Sx = Vx.*t_step;
    Sy = Vy.*t_step;
    
    x_pos = prev_x_pos + Sx;
    y_pos = prev_y_pos + Sy;

    
    for k = 1:display_particle
        subplot(211)
        plot([prev_x_pos(k) x_pos(k)],[prev_y_pos(k) y_pos(k)],'color',col(k,:));
  
        if x_pos(1,k) <= 0 
            x_pos(1,k) = boundary_x;
            y_pos(1,k) = prev_y_pos(1,k) ;
        elseif x_pos(1,k) >= boundary_x
            x_pos(1,k) = 0;
            y_pos(1,k) = prev_y_pos(1,k);
        end    
    
        if y_pos(1,k) <= 0 
            y_pos(1,k) = 0;
            Vy(1,k) = -(Vy(1,k));

        elseif  y_pos(1,k) >= boundary_y
            y_pos(1,k) = boundary_y;
            Vy(1,k) = -(Vy(1,k));
        end
        hold on

% Question 2
% Model the scattering of the electrons using an exponential scattering probability
        if Pscat > rand()
            Vx(1,k) = initial_Vth*randn;
            Vy(1,k) = initial_Vth*randn;
            VTot(1,k) = sqrt(Vx(1,k).^2 + Vy(1,k).^2);
        end
           
    end
    pause(0.001); 
    
    meanTemperature = mean(((VTot.^2).* mass_eff)./(2*C.kb));
    
    % Question 3
    % Average temperature over time
    subplot(212)
    plot([prev_time time], [prev_meanTemperature meanTemperature],'r')
    ylim([100 500]);
    xlabel('Time (s)');
    ylabel('Temperature (K)');
    hold on
    
    prev_meanTemperature = meanTemperature;
    
    pause(0.01)
    
    
end


