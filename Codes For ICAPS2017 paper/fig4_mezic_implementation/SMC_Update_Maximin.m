function [posagents, Ck] = SMC_Update_Maximin(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed)
    
    Lx = DomainBounds.xmax - DomainBounds.xmin;
    Ly = DomainBounds.ymax - DomainBounds.ymin;
    xmin = DomainBounds.xmin;
    ymin = DomainBounds.ymin;
    
    Nkx = size(muk, 1);
    Nky = size(muk, 2);
    Nagents = size(posagents, 1);
    
    persistent drag dist_thresh
    if isempty(dist_thresh)
        dist_thresh = 0.1;
    end
    if isempty(drag)
        drag = 0;
    end
    
    if drag>40
        drag = 0;
        dist_thresh = dist_thresh - 0.05; 
    end

    for iagent = 1:Nagents
        xrel = posagents(iagent, 1) - xmin;
        yrel = posagents(iagent, 2) - ymin;
        for kx = 0:Nkx-1
            for ky = 0:Nky-1        
    %             hk = sqrt(Lx*Ly); %hadi
                hk = Lx*Ly; %hadi
                if kx ~= 0
                    hk = hk * 0.5;
                end
                if ky ~= 0
                    hk = hk * 0.5;
                end
                hk = sqrt(hk);%hadi        
                Ck(kx+1, ky+1) = Ck(kx+1, ky+1) + (1/hk)*cos(kx * pi * xrel/Lx) * cos(ky * pi * yrel/Ly) * dt;
            end
        end
    end

    % Bhaskar - for all agents
    dist_agents = ones(Nagents, Nagents);
    for i = 1:Nagents
        for j = 1:Nagents
            if norm(posagents(i,:)-posagents(j,:))>dist_thresh
                dist_agents(i,j) = 0; 
            end
        end
    end 

    % Computing controls
    s=1.5;
    posagents_backup = posagents; 
    for iagent = 1:Nagents
        prtcpt_agents = posagents_backup(dist_agents(i,:)==1,:);
        if (length(prtcpt_agents(:,1))==1)        
            %Ergodicity_Metric=0;
            Bjx = 0.0; %hadi should be inside loop not outside
            Bjy = 0.0; %hadi
            
            xrel = posagents(iagent, 1) - xmin;
            yrel = posagents(iagent, 2) - ymin;
            for kx = 0:Nkx-1
                for ky = 0:Nky-1
                    lambda_k = 1.0 / ((1.0 + kx * kx + ky * ky)^s);
        %               hk = sqrt(Lx*Ly); %hadi
                      hk = Lx*Ly; %hadi
                    if kx ~= 0
                        hk = hk * 0.5;
                    end
                    if ky ~= 0
                        hk = hk * 0.5;
                    end
                    hk = sqrt(hk); %hadi
        
                    Bjx = Bjx + (lambda_k / hk) * (Ck(kx+1, ky+1) - Nagents*time*muk(kx+1, ky+1)) * (-kx * pi/Lx) * sin(kx * pi * xrel/Lx) * cos(ky * pi * yrel/Ly);
                    Bjy = Bjy + (lambda_k / hk) * (Ck(kx+1, ky+1) - Nagents*time*muk(kx+1, ky+1)) * (-ky * pi/Ly) * cos(kx * pi * xrel/Lx) * sin(ky * pi * yrel/Ly);
                end
            end
            
            Bjnorm = sqrt(Bjx*Bjx + Bjy*Bjy);
        
            % Updating agent location based on SMC feedback control law
            posagents(iagent,:) = posagents(iagent,:) - AgentSpeed * ([Bjx,Bjy]./Bjnorm) * dt;
            
            %No Need for this
            % reflecting agent in case it goes out of domain bounds
            posagents(iagent, :) = reflect_agent(posagents(iagent, :), DomainBounds);
        else
            if drag == 0
                dist_thresh = dist_thresh + 0.05; 
            end
            drag = drag + 1; 
            centroid = sum(prtcpt_agents)./length(prtcpt_agents(:,1));
            Bj = centroid - posagents(iagent,:);
            posagents(iagent,:) = posagents(iagent,:) - AgentSpeed * ([Bj(1),Bj(2)]./norm(Bj)) * dt;
            
            %No Need for this
            % reflecting agent in case it goes out of domain bounds
            posagents(iagent, :) = reflect_agent(posagents(iagent, :), DomainBounds);
        end
    end

%     % Bhaskar - only for two agents
%     dir_agents = posagents(1,:)-posagents(2,:); 
%     if norm(dir_agents)<0.3
%         posagents(1,:) = posagents(1,:) + AgentSpeed * ([dir_agents(1),dir_agents(2)]./norm(dir_agents)) * dt;
%         posagents(1, :) = reflect_agent(posagents(1, :), DomainBounds);
%         posagents(2,:) = posagents(2,:) - AgentSpeed * ([dir_agents(1),dir_agents(2)]./norm(dir_agents)) * dt;
%         posagents(2, :) = reflect_agent(posagents(2, :), DomainBounds);
%     else
%         % Updating Fourier Coefficients of Coverage Distribution
%         for iagent=1:Nagents
%             
%         end
%     end
end

function [agent] = reflect_agent(agent, DomainBounds)

    xmin = DomainBounds.xmin;
    xmax = DomainBounds.xmax;
    ymin = DomainBounds.ymin;
    ymax = DomainBounds.ymax;

    if agent(1) < xmin
        agent(1) = xmin + (xmin - agent(1, 1));
    end
    if agent(1) > xmax
        agent(1) = xmax - (agent(1, 1) - xmax);
    end
    if agent(2) < ymin
        agent(2) = ymin + (ymin - agent(1, 2));
    end
    if agent(2) > ymax
        agent(2) = ymax - (agent(1, 2) - ymax);
    end

end
