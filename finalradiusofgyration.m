function [n25,n10,n5]=finalradiusofgyration()
n25=radiusgyration5(25)
n10=radiusgyration5(10)
n5=radiusgyration5(5)
steps=200000

visualize_rg (n25,n10,n5,steps); %visualize conformation of the protein
end


function visualize_rg(n25,n10,n5,steps) %visualize N spheres at positions x

plot(n25(:,1),n25(:,2),n10(:,1),n10(:,2),n5(:,1),n5(:,2)); 
%xlim([0 10000])
%ylim([1.5 2.5])
%hold on; 
%view(2); hold on; 
end


function results=radiusgyration5(N)

close all;
set(0,'Defaultlinelinewidth',5, 'DefaultlineMarkerSize',6,...
    'DefaultTextFontSize',5, 'DefaultAxesFontSize',18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize parameters 

% # of amino acids
T=0.05; % temperature
dt=0.0005; % integration time step
steps=200000; % # of time steps

epsilon_LJ = 1; %LJ potential strength 
spring_coeff = 5; % Spring potential strength
sigma = 1.122; % equilibrium distance for both the LJ and spring potential 

L = sigma*N; % system size for visualization 
print_interval = 1000; % how often to spring the protein's new conformation  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize position velocities and forces  
x = initial_configuration(sigma,N); %init. coordinates x,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 

%results=zeros(400,2);

for step_i=1:steps %dynamics loop
    [x,pairs] = steepest_descent (N,x,dt,epsilon_LJ,sigma,spring_coeff,T);%calculate new positions based on forces 
    [ip,pair,connector]= all_interactions (N,x);
    out=radiusg(N,ip,connector);
    %answerA=step_i;
    %answerB=out;
    results(step_i,:)=[step_i out];
    %if mod(step_i-1,print_interval)==0
        %mytitle = ['step=',num2str(step_i), ' N=',num2str(N), ' L=',num2str(L)] ;
        %visualize_rg (results,steps); %visualize conformation of the protein
    %end
end
y_out = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rg=radiusg(N,ip,connector)
r=0;
for i=1:ip
    r=r+sqrt(((connector(i,1)).^2)+((connector(i,2)).^2));
    
end
rg=(1/(2.*(N.^2))).*r
end


function x= initial_configuration (initial_sigma,N)
x = zeros(N,2); 
x(:,1) = initial_sigma.*[1:1:N] - (initial_sigma*N/2); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,pairs] = steepest_descent (N,x,dt,epsilon_LJ,sigma,spring_coeff,T) 
[F_particles,~,pairs] = forces(N,x,epsilon_LJ,sigma,spring_coeff);
F = F_particles; 
x = x + (dt*F)+ T.*(rand(size(x))-0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= all_interactions (N,x) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,2);
for i=1:N-1
    for j=i+1:N
        distance = (x(j,:)-x(i,:));
        ip = ip + 1; %interaction pair counter
        pair(ip,:) = [i j]; %particle numbers (i,j) belonging to pair (ip)
        connector(ip,:) = distance;
    end
end %end both for loops
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= spring_interactions (N,x) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,2); 
for i=1:N-1
    j=i+1;
    distance = (x(j,:)-x(i,:));
    ip = ip + 1; %interaction pair counter
    pair(ip,:) = [i j]; %particle numbers (i,j) belonging to pair (ip)
    connector(ip,:) = distance;
end %end for loops
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,P,pair] = forces (N,x,epsilon_LJ,sigma,spring_coeff) %clear forces F, then calculate them using ..
F=zeros(N,2); P=zeros(N,2);
% LJ forces 
[no,pair,connector]=all_interactions(N,x); %interacting pairs
for i=1:no, FORCE=force_LJ(connector(i,:),epsilon_LJ);
    F(pair(i,1),:)=F(pair(i,1),:)-FORCE; F(pair(i,2),:)=F(pair(i,2),:)+FORCE; %actio=reactio; 
    P(pair(i,1),:)=P(pair(i,1),:)+(sum(FORCE.*connector(i,:))); 
    P(pair(i,2),:)=P(pair(i,2),:)+(sum(FORCE.*connector(i,:)));
end
%spring forces:
[no,pair,connector]=spring_interactions(N,x); %interacting pairs
for i=1:no
    FORCE = force_springs(connector(i,:),spring_coeff,sigma);
    F(pair(i,1),:)=F(pair(i,1),:)-FORCE; F(pair(i,2),:)=F(pair(i,2),:)+FORCE; %actio=reactio; 
    P(pair(i,1),:)=P(pair(i,1),:)+(sum(FORCE.*connector(i,:))); 
    P(pair(i,2),:)=P(pair(i,2),:)+(sum(FORCE.*connector(i,:)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function curr_force = force_springs(r_vector,spring_coeff_array,sigma)
r2=sum(r_vector.^2,2);
r = sqrt(r2);
curr_force = zeros(length(r2), 2);
curr_force(:,1) = -spring_coeff_array'.*((r -  sigma)).*(r_vector(:,1)./r);
curr_force(:,2) = -spring_coeff_array'.*((r -  sigma)).*(r_vector(:,2)./r);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force_LJ= force_LJ (r_vector,epsilon_LJ); r=norm(r_vector); 
force_LJ = 24*epsilon_LJ*(2*r.^(-14)-r^(-8)) * r_vector; %Lennard Jones
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visualize_particles (N,x,L,pairs,mytitle) %visualize N spheres at positions x

%plot particles:
for i=1:N
    scatter(x(i,1),x(i,2),80,'b','filled'); hold on; 
    view(2); hold on; 
end

%plot springs:
for index_i=1:length(pairs)
    curr_index =  index_i;
    curr_i = pairs(curr_index,1);
    curr_j = pairs(curr_index,2);
    
    plot([x(curr_i,1),x(curr_j,1)], ...
        [x(curr_i,2),x(curr_j,2)],...
        'Color',[.6,.6,.6],'LineWidth',1.0); hold on;
end

%colorbar; 
xlim([-L,L]/2); ylim([-L,L]/2);
axis square; hold off; 
title(mytitle); pause(0.001); 
end



%colorbar; 
%xlim([-L,L]/2); ylim([-L,L]/2);
%axis square; hold off; 
%title(mytitle); pause(0.001); 
%end