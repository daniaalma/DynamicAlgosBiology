function out=random_walks()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feb 2019, Orit Peleg, orit.peleg@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
set(0,'Defaultlinelinewidth',2, 'DefaultlineMarkerSize',3,...
    'DefaultTextFontSize',15, 'DefaultAxesFontSize',15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=500; % number of steps per trajectory 
REALIZATIONS = 50; %number of trajectories
v = 1; % velocity (aka step size)
theta_s_array =[pi/24, pi/12, pi/3]; %[pi/24,pi/6,pi/3]; % the width of the random walk turning angle distribution (the lower it is, the more straight the trajectory will be)  
w_array =0:0.01:1;%[0.0,0.5,1]; %0:0.01:1; % w is the weighting given to the directional bias (and hence (1-w) is the weighting given to correlated motion) 
ratio_theta_s_BRW_CRW = 0.5;  
PLOT_WALKS=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

efficiency_array = zeros(length(theta_s_array),length(w_array));

for w_i=1:length(w_array)
    w = w_array(w_i);
    for theta_s_i = 1:length(theta_s_array)
        theta_s_CRW = ratio_theta_s_BRW_CRW*theta_s_array(theta_s_i);
        theta_s_BRW = theta_s_array(theta_s_i);
        [x,y] = BCRW(N,REALIZATIONS, v, theta_s_CRW,theta_s_BRW,w);
        if PLOT_WALKS==1
            figure();
            plot(x',y');axis equal;  xlim([0 N]);  hold on;
        end
        efficiency_array(theta_s_i,w_i) = mean(x(:,end)-x(:,1))./(v*N); 
    end
end

figure(); legend_array = []; 
for theta_s_i = 1:length(theta_s_array)
    theta_s_CRW = ratio_theta_s_BRW_CRW*theta_s_array(theta_s_i);
    theta_s_BRW = theta_s_array(theta_s_i);
    legend_array{theta_s_i} = ...
        ['\theta^{*CRW}=',num2str(theta_s_CRW,'%.2f'),...
        ' \theta^{*BRW}=',num2str(theta_s_BRW,'%.2f')]; 
    plot (w_array,efficiency_array(theta_s_i,:),'-o'); hold on; 
    
end
hold off; legend(legend_array); xlabel('w'); ylabel('navigational efficiency'); 

out = 0;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this funciton generate 2D Biased Corrolated Random Walks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=BCRW(N,REALIZATIONS,v,theta_s_CRW,theta_s_BRW,w)

x = zeros(REALIZATIONS, N);
y = zeros(REALIZATIONS, N);
theta = zeros(REALIZATIONS, N);
x(:,1) = 0;
y(:,1) = 0;
theta(:,1) = 0;

for REALIZATION_i = 1:REALIZATIONS
    for step_i = 2:N
        theta_CRW = theta(REALIZATION_i,step_i - 1)+(theta_s_CRW*2.0*(rand(1,1)-0.5));
        theta_BRW = theta_s_BRW*2.0*(rand(1,1)-0.5);
 
        x(REALIZATION_i,step_i) = x(REALIZATION_i,step_i -1) + ...
            v*((w*cos(theta_BRW)) + ((1-w)*cos(theta_CRW)));
        y(REALIZATION_i,step_i) = y(REALIZATION_i,step_i -1) + ...
            v*((w*sin(theta_BRW)) + ((1-w)*sin(theta_CRW)));
        
        current_x_disp = x(REALIZATION_i,step_i)-x(REALIZATION_i,step_i-1);
        current_y_disp = y(REALIZATION_i,step_i)-y(REALIZATION_i,step_i-1); 
        current_direction = atan2(current_y_disp,current_x_disp); 
        
        theta(REALIZATION_i,step_i) = current_direction;
        
    end
end
end








