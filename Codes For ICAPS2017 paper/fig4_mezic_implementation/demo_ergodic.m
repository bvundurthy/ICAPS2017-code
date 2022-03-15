close all
clear variables
clc

%% Setting domain bounds
DomainBounds.xmin = 0.0;
DomainBounds.xmax = 1.0;
DomainBounds.ymin = 0.0;
DomainBounds.ymax = 1.0;
Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;

%obstacles
obstacles.r = [];%[0.1,0.05,0.08,0.05];% radii of obstacles
obstacles.p = [];%[[0.5;0.5],[0.2;0.3],[0.7;0.2],[0.3;0.8]];%positions of obstacles[x1 y1
%                                                                               x2 y2
                                                                               .....]
obstacles.number = numel(obstacles.r);

% Number of wave-numbers to be used
Nk = 50;%%
livePlot = true; %set <true> if you want live plot for trajectories, other <false> for faster execution
%% Calculating muk
res = 100;% resolution of discretization
mu=ones(res,res);
xdel=Lx/res;
ydel=Ly/res;


% for xRange=0:xdel:Lx-xdel
%          for yRange=0:ydel:Ly-ydel
%             for i=1:obstacles.number
%                 if (xRange-obstacles.p(1,i))^2 + (yRange-obstacles.p(2,i))^2 <= obstacles.r(i)^2
%                     mu(uint8(xRange*res+1),uint8(yRange*res+1)) = 0;
%                 end
%             end
%          end
% end
mu=mu./sum(sum(mu));
[X,Y]=meshgrid(1:res,1:res);
surf(X,Y,mu);
close;
% muk = dct2(mu,Nk,Nk);
%For Matlab DCT to work: (k1,k2)=1,2,..(L1-L2)

muk=zeros(Nk,Nk);
xdel=Lx/res;
ydel=Ly/res;
Nkx = size(muk, 1);
Nky = size(muk, 2);
for kx = 0:Nkx-1
    for ky = 0:Nky-1
        
        hk=Lx*Ly; %using lim x->0 sinx/x=1
        if kx ~= 0
            hk = hk * 0.5;
        end
        if ky ~= 0
            hk = hk * 0.5;
        end
        hk = sqrt(hk);
        
        for xRange=0:xdel:Lx-xdel
            for yRange=0:ydel:Ly-ydel
                muk(kx+1, ky+1) = muk(kx+1, ky+1)+ mu(uint8(xRange*res+1),uint8(yRange*res+1)) *(1/hk)*cos(kx * pi * xRange/Lx) * cos(ky * pi * yRange/Ly);
                
            end
        end
        
        
    end
end

%% Initializing agent locations
Nagents = 1;%4;
posagents = [0.35,0.7];%[0.35,0.7;0.05,0.3;0.1,0.95;0.83,0.03];%make sure not to start in an obstacle region!
AgentSpeed = 25;

colors = ['m','g','b','c'];

Nsteps = 2500;
dt = 0.001;

% Initializing Fourier coefficients of coverage distribution
Ck = zeros(Nk, Nk);
ck_t = zeros(Nk, Nk);

viscircles(obstacles.p',obstacles.r);
%color inside the circles
for xRange=0:xdel:Lx-xdel
         for yRange=0:ydel:Ly-ydel
            for i=1:obstacles.number
                if (xRange-obstacles.p(1,i))^2 + (yRange-obstacles.p(2,i))^2 <= obstacles.r(i)^2
                    scatter(xRange, yRange,2,'r','fill');
                end
            end
         end
end

axis equal;

phi_squared = zeros(Nsteps,1);
s= 1.5;
% Executing multiple steps of SMC algorithm
Ergodicity_Metric_save=0;
tic

traj_agent1 = zeros(Nsteps,2); 
for it = 1:Nsteps
    time = (it) * dt;
    [posagents, Ck] = SMC_Update(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed);
    ck_t = Ck/(Nagents*time);    
    traj_agent1(it,:) = posagents(1, :);
%     for iagent = 1:Nagents
%         traj_agent1 = posagents(iagent, :);
%         if livePlot == true
%                 pause(0.0001)
%         end
%     end
%     if mod(it,50)==0     
%         fprintf('Iteration: %i/%i  \n', it,Nsteps) 
%     end
    
    [Ergodicity_Metric] = Calculate_Ergodicity(Ck/Nagents/time, muk,DomainBounds);
    Ergodicity_Metric_save=[Ergodicity_Metric_save,Ergodicity_Metric];
    
end
toc


%% plotting part
close all; 
clc;

figure(1); 
subplot(1,2,1); axis ([-0.1 1.1 -0.1 1.1]); hold on;
p2 = subplot (1,2,2); hold on; xlim([0 100]); ylim([0 100]); zlim([0 25]); 

cove = zeros(res, res); 

for it = 2:Nsteps
    subplot(1,2,1);
    line([traj_agent1(it-1, 1) traj_agent1(it, 1)], [traj_agent1(it-1, 2) traj_agent1(it, 2)], 'Color', colors(1), 'Marker', 'o', 'MarkerSize', 1);        

    iagent_x = round(traj_agent1(it, 1)*res); iagent_y = round(traj_agent1(it, 2)*100);
    iagent_xmin = max(iagent_x-3,1); iagent_xmax = min(iagent_x+3,res);
    iagent_ymin = max(iagent_y-3,1); iagent_ymax = min(iagent_y+3,res);
    cove(iagent_xmin:iagent_xmax, iagent_ymin:iagent_ymax) = cove(iagent_xmin:iagent_xmax, iagent_ymin:iagent_ymax) +1; 

    pause(0.00001);     
    p2 = subplot(1,2,2);
    surf(Y,X,cove);
    if mod(it,50)==0 && it~=Nsteps
        cla(p2); 
        fprintf('Iteration: %i/%i  \n', it,Nsteps) 
    end
end
% 
% figure
% surf(X,Y,cove); 

%% Plotting the metric of ergodicity 
% time=0:0.001:0.001*Nsteps;
% figure;loglog(time(2:end),Ergodicity_Metric_save(2:end))
% axis([0.001 5 0.0001,1])
% xlabel('Time');
% ylabel('Coverage Metric, \phi(t)');
% title('Metric of ergodicity as a function of time')


figure

