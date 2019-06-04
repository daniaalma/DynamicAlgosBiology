
function y_out = Flocking_MatlabF()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MArch 2019, Orit Peleg, orit.peleg@colorado.edu
% Code for HW3 CSCI 4314/5314 Dynamic Models in Biology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
set(0,'Defaultlinelinewidth',5, 'DefaultlineMarkerSize',6,...
    'DefaultTextFontSize',5, 'DefaultAxesFontSize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize parameters
N = 200;        %No. of boids
frames = 100;   %No. of frames in movie
limit = 50;    %Axis limits
L=limit*2; 
P = 10;         %Spread of initial position (gaussian)
V = 10;         %Spread of initial velocity (gaussian)
delta = 1; %0.01;  %Time step
c1 = 0.00001;       %Attraction scaling factor
c2 = 0.001;         %Repulsion scaling factor
c3 = 1;          %Heading scaling factor
c4 = 0.01;         %Randomness scaling factor
vlimit = 1;    %Maximum velocity


%Initialize
p = P*randn(2,N);
v = V*randn(2,N);
figure(); 


%Main loop
for k=1:frames
    v1=zeros(2,N);
    v2=zeros(2,N);
    v3 = [sum(v(1,:))/N; sum(v(2,:))/N]*c3;             %Calculate average veolcity
    if(norm(v3) > vlimit), v3 = v3*vlimit/norm(v3); end %Limit max velocity
    for n=2:N
        for m=2:N %%%
            if m ~= n
                r = p(:,m)-p(:,n);                  %Vector from one agent to the next
                
                dist=p(:,m)-p(:,1);
                
                if r(1)>L/2, r(1) = r(1)-L;
                elseif r(1)<-L/2, r(1) = r(1)+L;
                end
                
                if r(2)>L/2, r(2) = r(2)-L;
                elseif r(2)<-L/2, r(2) = r(2)+L;
                end
                
                rmag=sqrt(r(1)^2+r(2)^2);           %Distance between agents
                v1(:,n) = v1(:,n) + c1*r;           %Attraction
                v2(:,n) = v2(:,n) - c2*r/(rmag^2);  %Repulsion (non-linear scaling)
                
                if dist<10 & dist>-10, %dist=dist+20;%%%

                    v2(:,n) = v2(:,n) - 300*r/(rmag^2);
                end
                
            end
        end
        v4(:,n) = c4*randn(2,1);                    %Add randomness to motion
        v(:,n) = v1(:,n)+v2(:,n)+v4(:,n)+v3;        %Update velocity
    end
    p = p + v *delta;            %Update position
    p(:,1)=[0 0];
    %%%
    
    % Periodic boundary 
    tmp_p = p;
    
    
    tmp_p(1,p(1,:)>L/2) = tmp_p(1,p(1,:)>L/2) - L;
    tmp_p(2,p(2,:)>L/2) = tmp_p(2,p(2,:)>L/2) - L;
    tmp_p(1,p(1,:)<-L/2)  = tmp_p(1,p(1,:)<-L/2) + L;
    tmp_p(2,p(2,:)<-L/2)  = tmp_p(2,p(2,:)<-L/2) + L;
    
    p = tmp_p;
    
    %Update plot:
    plot(p(1,:),p(2,:),'k+','Markersize',4); 
    quiver(p(1,:),p(2,:),v(1,:),v(2,:)); %For drawing velocity arrows
    axis([-limit limit -limit limit]); axis square; 
    drawnow; 
    hold off;
end
y_out=0; 


end
