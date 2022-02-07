%% Part 4 Injection
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
C.g = 9.80665;                      % metres (32.1740 ft) per sÂ²

boundary_x = 200e-9;                  
boundary_y = 100e-9;
t_step = 5e-15;
tao_mn = 0.2e-12;
n_particle = 1e4;

turnOnScatter = 1;                  % Switch on to turn on electron scattering

inject_particle = 2;                % Number of inject particles
injection_Period = 20;              % Injection time period
injection_x =  0.01e-7;             % Electron injection x position
injection_y_max = 0.52e-7;          % Electron injection y position
injection_y_min = 0.48e-7;          % Electron injection y position
nSims = 500;
T = 300;

mass_eff =  0.26*C.m_0;
initial_Vth = sqrt((2* C.kb *T)/mass_eff);



% Generate electrons x and y velocity
Vx = initial_Vth.*abs(randn(1,inject_particle));        
% Vy = initial_Vth.*randn(1,inject_particle);
Vy = (Vx.*0.5).*rand(1,inject_particle);

% Generate electron with the same x positions
x_pos = zeros(1,inject_particle)+injection_x;   
% Generate electron y positions within a certian range
y_pos = (injection_y_max - injection_y_min).*rand(1,inject_particle)+injection_y_min;


Pscat = 1- exp(-(t_step/tao_mn));

display_particle = inject_particle; % display the number of injected electrons

% Generate the space of 200e-7 * 100e-7 and three boxes

figure(2)
rectangle('position',[0.8e-7 0e-7 0.4e-7 0.4e-7]);
rectangle('position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7]);
rectangle('position',[1.4e-7 0.4e-7 0.2e-7 0.2e-7]);
hold on
%col = hsv(10);
xlim([0 200e-9]);
ylim([0 100e-9]);





for n = 1: nSims
    %% Introduce them during the simulation from the left side with a positive vx –
    %  derived from a thermalized velocity within a small central region.
    
    
    % Every injection period, introduce 2 electrons into the area
    if n >= injection_Period && mod(n,injection_Period) == 0
        new_inject_x = zeros(1,inject_particle)+0.01e-7;
        new_inject_y = (injection_y_max - injection_y_min).*rand(1,inject_particle)+injection_y_min;
        x_pos = [x_pos new_inject_x];
        y_pos = [y_pos new_inject_y];
        % Generate new electrons' x and y positions and add them to the
        % x and y position array
        
        new_inject_Vx = initial_Vth.*abs(randn(1,inject_particle));
        new_inject_Vy = initial_Vth.*randn(1,inject_particle)-0.5;
        Vx = [Vx new_inject_Vx];
        Vy = [Vy new_inject_Vy];
        display_particle = display_particle + inject_particle;
        % Generate new electrons velocity on x and y direction and add them to the
        % x and y velocity arrays
    end
    
    prev_x_pos = x_pos;         % save the old x position       
    prev_y_pos = y_pos;         % save the old y position
    
    Sx = Vx.*t_step;            % Update x position = X velocity * time step 
    Sy = Vy.*t_step;            % Update y position = Y velocity * time step 
    
    x_pos = prev_x_pos + Sx;    % New x position 
    y_pos = prev_y_pos + Sy;    % New y position 

    col = hsv(display_particle);     % Color of the electron movement
    for k = 1:display_particle
        plot([prev_x_pos(1,k) x_pos(1,k)],[prev_y_pos(1,k) y_pos(1,k)],'color',col(k,:));
        % plot all electrons within the same time step 
 
        %% Turn off the periodic BC conditions in x
   
        % Left X boundary condition
        if x_pos(1,k) <= 0 
            x_pos(1,k) = 0;
            prev_x_pos(1,k) = 0;
            Vx(1,k) = -(Vx(1,k));
            
        % Right X boundary condition
        elseif x_pos(1,k) >= boundary_x
            Vx(1,k) = -(Vx(1,k));
            x_pos(1,k) = boundary_x;
            prev_x_pos(1,k) = boundary_x;
            
        % Bottom Y boundary condition    
        elseif y_pos(1,k) <= 0 
            y_pos(1,k) = 0;
            Vy(1,k) = -(Vy(1,k));
            prev_y_pos(1,k) = 0;

        % Top Y boundary condition  
        elseif  y_pos(1,k) >= boundary_y
            y_pos(1,k) = boundary_y;
            prev_y_pos(1,k) = boundary_y;
            Vy(1,k) = -(Vy(1,k));
        
        % 1st and 2nd box left side boundary condition  
        elseif (x_pos(1,k) >= 0.8e-7 && x_pos(1,k) <= 1e-7 ) && (y_pos(1,k) <= 0.4e-7 || y_pos(1,k) >= 0.6e-7) && prev_x_pos(1,k) < 0.8e-7
            x_pos(1,k) = 0.8e-7;
            prev_x_pos(1,k) = 0.8e-7;
            Vx(1,k) = -(Vx(1,k));
            
        % 1st and 2nd box right side boundary condition  
        elseif (x_pos(1,k) >= 1e-7 && x_pos(1,k) <= 1.2e-7) && (y_pos(1,k) <= 0.4e-7 || y_pos(1,k) >= 0.6e-7) && prev_x_pos(1,k) <= 1.2e-7  
            x_pos(1,k) = 1.2e-7;
            prev_x_pos(1,k) = 1.2e-7;
            Vx(1,k) = -(Vx(1,k));


        % 1st box top side boundary condition  
        elseif (y_pos(1,k) <= 0.4e-7) && (x_pos(1,k) >= 0.8e-7 && x_pos(1,k) <= 1.2e-7)
            y_pos(1,k) = 0.4e-7;
            prev_y_pos(1,k) = 0.4e-7;
            Vy(1,k) = -(Vy(1,k));
            
        % 2nd box bottom side boundary condition  
        elseif (y_pos(1,k) >= 0.6e-7) && (x_pos(1,k) >= 0.8e-7 && x_pos(1,k) <= 1.2e-7)
            y_pos(1,k) = 0.6e-7;
            prev_y_pos(1,k) = 0.6e-7;
            Vy(1,k) = -(Vy(1,k));
%% Adding the third new box 
            
%       3rd box left side boundary condition  
        elseif (y_pos(1,k) >= 0.4e-7 && y_pos(1,k) <= 0.6) && (x_pos(1,k) >=1.4e-7 && x_pos(1,k) <= 1.5e-7) && prev_x_pos(1,k)<1.4e-7
            x_pos(1,k) = 1.4e-7;
            prev_y_pos(1,k) = 1.4e-7;
            Vx(1,k) = -(Vx(1,k));   

%       3rd box right side boundary condition  
        elseif (y_pos(1,k) >= 0.4e-7 && y_pos(1,k) <= 0.6) && (x_pos(1,k) >=1.5e-7 && x_pos(1,k) <= 1.6e-7)&& prev_x_pos(1,k)>1.6e-7
            x_pos(1,k) = 1.6e-7;
            prev_y_pos(1,k) = 1.6e-7;
            Vx(1,k) = -(Vx(1,k));   

%       3rd box top side boundary condition  
        elseif (y_pos(1,k) <= 0.6e-7 && y_pos(1,k) >= 0.5e-7) && (x_pos(1,k) >=1.4e-7 && x_pos(1,k) <= 1.6e-7) && prev_y_pos(1,k)>0.6e-7
            y_pos(1,k) = 0.6e-7;
            prev_y_pos(1,k) = 0.6e-7;
            Vy(1,k) = -(Vy(1,k));   

%       3rd box bottom side boundary condition  
        elseif (y_pos(1,k) >= 0.4e-7 && y_pos(1,k) <= 0.5e-7) && (x_pos(1,k) >=1.4e-7 && x_pos(1,k) <= 1.6e-7) && prev_y_pos(1,k)<0.4e-7
            y_pos(1,k) = 0.4e-7;
            prev_y_pos(1,k) = 0.4e-7;
            Vy(1,k) = -(Vy(1,k));   
        end 
        
        hold on
        
        % if the generated random number bigger than the scattering
        % probability number, the electron scatters
        if turnOnScatter == 1           
            if Pscat > rand()
                Vx(1,k) = initial_Vth*randn;
                Vy(1,k) = initial_Vth*randn;
            end
        end
        
          
    end
    pause(0.001); 
    
end
