%% Part 3 Enhancement
% Qiushi Chen 101049864


clc
clear all
close all

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
t_step = 5e-15;

n_particle = 1e4;
display_particle = 10;
nSims = 500;
T = 300;
loop_index = 10000;

% Question 2
% Make the boundary diffusive by setting a velocity loss parameter
Vloss = 0.8;                % Velocity loss parameter
t_mn = 0.2e-12;
mass_eff =  0.26*C.m_0;
initial_Vth = sqrt((2* C.kb *T)/mass_eff);

x_pos = zeros(1,display_particle);
y_pos = zeros(1,display_particle);

n_x_pos = boundary_x*rand(1,n_particle);
n_y_pos = boundary_y*rand(1,n_particle);

% Generate the electron positions outside two boxes
index = 1;
for k = 1:loop_index
    rand_x = 2*rand;
    rand_y = rand;
    
    % Check whether the x and y positions are outside the box
    if(rand_x < 0.8 || rand_x > 1.2) || ((rand_y > 0.4 && rand_y < 0.6) && (rand_x >= 0.8 || rand_x <= 1.2))
        x_pos(1,index) = 1e-7*rand_x; 
        y_pos(1,index) = boundary_y*rand_y;
        index = index + 1;
    end
    if index == display_particle + 1
        break % Once the number of electrons satisfies the requirement, break the loop
    end
end    


Vx = initial_Vth.*randn(1,display_particle);
Vy = initial_Vth.*randn(1,display_particle);

Pscat = 1- exp(-(t_step/t_mn));

% Question 1
% Add in the inner rectangle "bottle neck" boundaries
figure(2)
rectangle('position',[0.8e-7 0e-7 0.4e-7 0.4e-7]);
rectangle('position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7]);
hold on

col = hsv(10);          % Color of trajectory
xlim([0 200e-9]);
ylim([0 100e-9]);

xlabel('Length (m)');
ylabel('Width (m)');


for n = 1: nSims
    
    prev_x_pos = x_pos;
    prev_y_pos = y_pos;
   
    
    
    Sx = Vx.*t_step;
    Sy = Vy.*t_step;
    
    x_pos = prev_x_pos + Sx;
    y_pos = prev_y_pos + Sy;

    
    for k = 1:display_particle
        %subplot(211)
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
            Vy(1,k) = -(Vy(1,k))*Vloss;     % Boundary diffusive

        elseif  y_pos(1,k) >= boundary_y
            y_pos(1,k) = boundary_y;
            Vy(1,k) = -(Vy(1,k))*Vloss; % Boundary diffusive
        end
        
        if (x_pos(1,k) >= 0.8e-7 && x_pos(1,k) <= 1e-7 ) && (y_pos(1,k) <= 0.4e-7 || y_pos(1,k) >= 0.6e-7) && prev_x_pos(1,k) < 0.8e-7
            x_pos(1,k) = 0.8e-7;
            prev_x_pos(1,k) = 0.8e-7;
            Vx(1,k) = -(Vx(1,k));

        elseif (x_pos(1,k) >= 1e-7 && x_pos(1,k) <= 1.2e-7) && (y_pos(1,k) <= 0.4e-7 || y_pos(1,k) >= 0.6e-7) && prev_x_pos(1,k) > 1.2e-7  
            x_pos(1,k) = 1.2e-7;
            prev_x_pos(1,k) = 1.2e-7;
            Vx(1,k) = -(Vx(1,k));




        elseif (y_pos(1,k) <= 0.4e-7) && (x_pos(1,k) >= 0.8e-7 && x_pos(1,k) <= 1.2e-7)
            y_pos(1,k) = 0.4e-7;
            prev_y_pos(1,k) = 0.4e-7;
            Vy(1,k) = -(Vy(1,k))*Vloss;  % Boundary diffusive
        elseif (y_pos(1,k) >= 0.6e-7) && (x_pos(1,k) >= 0.8e-7 && x_pos(1,k) <= 1.2e-7)
            y_pos(1,k) = 0.6e-7;
            prev_y_pos(1,k) = 0.6e-7;
            Vy(1,k) = -(Vy(1,k))*Vloss;     % Boundary diffusive


         end 

        
        
        hold on
        
        if Pscat > rand()
            Vx(1,k) = initial_Vth*randn;
            Vy(1,k) = initial_Vth*randn;
        end
          
    end
    pause(0.001); 
    

    
   
  
    
end

% Question 3
% Electron density map
figure(3)
electronMapX = linspace(0,boundary_x, 100);
electronMapY = linspace(0,boundary_y, 50);
ElectonDensityMap = histcounts2(y_pos,x_pos,electronMapY,electronMapX);
%subplot(212);
imagesc(electronMapY,electronMapX,ElectonDensityMap),colorbar,title('Electron density map');
xlabel('Length (m)');
ylabel('Width (m)');



