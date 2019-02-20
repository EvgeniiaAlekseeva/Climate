% Dimensionality of the phenotypic space
D = 2;  

% Number of clusters
m = 1;  
m_plot = m;
mD = m*D;
m_max = 49;
% Phenotypes of initial clusters
X_end = normrnd(0,1,1,mD); 
X_plot = X_end; 
% Populations of initial clusters
N_end = 0.3.*ones(1,m);                               
% Location of CCC
x_cc = 0.*ones(1,D);  
% Width of capacity kernel
sigma_b = 1/sqrt(D);
% Width of competition kernel
sigma_comp = 0.5.*ones(1,D);                          
% Mutation coefficients
b = normrnd(0,sigma_b,D,D);                  

% Parameters of clusters evolution
dist_to_join = 0.001;   % critical distance to join clusters
dist_to_sep = 0.001;    % critical distance to split clusters
dist_to_plot = 0.1;     % criticial distance from which clusters are shown as separate on the plot
min_size = 0.000001;    % critical population size before disappearence

% Time perionds of ecological and evolutionaru dynamics
tspan_ecology = [0 500];
tspan_evolution = [0 0.1];
P = 700; % nuber of itterations
p = 1;
  
 
%% Plotting parameters and plot of initial conditions
% scale of axis on the plot
ax = 2;                                                
axx1= x_cc(1) - 2;
axx2 = x_cc(1) + 2;                                    
axy1 = x_cc(2) - 2;
axy2 = x_cc(2) + 2;
for k = 1:length(N_end)
        g = 2*k;
        size = 16*N_end(k)+4;
        set(0,'DefaultAxesFontSize',16);
        plot (X_end(g-1), X_end(g), 'or', 'MarkerSize', size, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'w');
        axis([axx1 axx2 axy1 axy2]);
        grid on;
        hold on;
end
plot (x_cc(1,1), x_cc(1,2), 'or', 'MarkerSize', 6, 'Marker','o','Color','r', 'LineWidth', 1.5);
set(gcf,'color','w');            
drawnow;                         
F(1) = getframe(gcf);            
pause(.1)


%% TIME STEPS CYCLE
while p < P
    hold off;
     
    %% CompetitionMatrix(i,j) counts how i-th cluster influense on j-th
    CompetitionMatrix = ones(m,m);
    for i = 1:m
        i1 = D*(i-1)+1;
        iD = i*D;
        for j = 1:m
            j1 = D*(j-1) +1;
            jD = j*D;
            CompetitionMatrix(i,j) = competitionD_cl(X_end(i1:iD), X_end(j1:jD), D, b, sigma_comp, x_cc);
        end
    end
    
    %% ECOLOGICAL STEP
    [t N] = ode45(@PopulDens_cl, tspan_ecology, N_end, [], D, b, m, X_end, CompetitionMatrix, x_cc);
    N_end = N(end, :);
    N = [];
    
    %% KILLING SMALL CLUSTERS
    for q = 1:length(N_end)
        if N_end(q) < min_size           
            N_end(q) = 0;
            X_end( (D*(q-1)+1):q*D ) = 0;
        end
    end
    N_end(N_end == 0) = [];              
    X_end(X_end == 0) = [];
    m = length(N_end);
    mD = m*D;
    X_end_previous = X_end;              
    
    %% EVOLUTIONARY STEP
    [t X] = ode45(@Phenotypes_cl, tspan_evolution, X_end, [], D, b, m, mD, sigma_comp, N_end, x_cc);
    X_end = X(end, :);
    X = [];
    
    %% JOINING CLOSE CLUSTERS
    if rem(p,10) == 5     
        cl = 1;
        while cl <= m
            clD1 = (cl-1)*D + 1;
            clD = cl*D;
            N_cl = N_end(cl);
            X_cl = X_end(clD1:clD);
            cl1 = cl + 1;
            while cl1 <= m
                cl1D1 = (cl1-1)*D + 1;
                cl1D = cl1*D;
                X_cl1 = X_end(cl1D1:cl1D);
                dist = sqrt(sum((X_cl1-X_cl).^2));
                if dist < dist_to_join
                    for dd = 1:D
                        X_cl(dd) = (X_cl(dd)*N_cl + X_cl1(dd)*N_end(cl1))/(N_cl+N_end(cl1));
                    end
                    N_cl = N_cl + N_end(cl1);
                    X_end(clD1:clD) = X_cl;
                    X_end(cl1D1:cl1D) = [];
                    N_end(cl) = N_cl;
                    N_end(cl1) = [];
                    m = m-1;
                    mD = m*D;
                else
                    cl1 = cl1 + 1;
                end
            end
            cl = cl +1;
        end
    end
     
    %% COUNTING VISIBLE CLUSTERS
    X_plot = X_end;
    N_plot = N_end;
    m_plot = m;
    pcl = 1;
        while pcl <= m_plot
            pclD1 = (pcl-1)*D + 1;
            pclD = pcl*D;
            N_pcl = N_plot(pcl);
            X_pcl = X_plot(pclD1:pclD);
            pcl1 = pcl + 1;
            while pcl1 <= m_plot
                pcl1D1 = (pcl1-1)*D + 1;
                pcl1D = pcl1*D;
                X_pcl1 = X_plot(pcl1D1:pcl1D);
                pdist = sqrt(sum((X_pcl1-X_pcl).^2));
                if pdist < dist_to_plot
                    for dd = 1:D
                        X_pcl(dd) = (X_pcl(dd)*N_pcl + X_pcl1(dd)*N_plot(pcl1))/(N_pcl+N_plot(pcl1));
                    end
                    N_pcl = N_pcl + N_plot(pcl1);
                    X_plot(pclD1:pclD) = X_pcl;
                    X_plot(pcl1D1:pcl1D) = [];
                    N_plot(pcl) = N_pcl;
                    N_plot(pcl1) = [];
                    m_plot = m_plot -1;
                else
                    pcl1 = pcl1 + 1;
                end
            end
            pcl = pcl +1;
        end
        
    
    %% PLOTTING OF THIS CYCLE'S FIGURE
    for k = 1:(length(N_plot))
        if k ==1
            g = 2*k;
            size = 16*N_plot(k)+4;
            plot (X_plot(g-1), X_plot(g), 'or', 'MarkerSize', size, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'w');
            set(0,'DefaultAxesFontSize',16);
            axis([axx1 axx2 axy1 axy2]);
            grid on;
            hold on;
        else
            g = 2*k;
            size = 16*N_plot(k)+4;
            plot (X_plot(g-1), X_plot(g), 'or', 'MarkerSize', size, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'w');
            set(0,'DefaultAxesFontSize',16);
            axis([axx1 axx2 axy1 axy2]);
            grid on;
            hold on;
        end
    end
    plot (x_cc(1,1), x_cc(1,2), 'or', 'MarkerSize', 6, 'Marker','o', 'Color', 'r', 'LineWidth', 1.5);
    xlabel('Phenotype 1','FontSize', 14, 'FontName', 'Serif', 'FontWeight', 'bold');
    ylabel('Phenotype 2','FontSize', 14, 'FontName', 'Serif', 'FontWeight', 'bold');
    title('Equillibration in D=2', 'FontName', 'Serif', 'FontWeight', 'normal');
    set(gcf,'color','w');           
    drawnow;                         
    if p > 0
        F(p) = getframe(gcf);            
        pause(.1)
    end
    
    
    %% MUTANT APPEARENCE
    if m < 7^D 
      mutant = randi ([1, m], 1, 1);           
      coef = rand(1,D);                        
      coef = coef./(sum(coef)/D);              
      eq_delta_x = dist_to_sep/(2*sqrt(D));
      delta_x = eq_delta_x.*sqrt(coef);   
      A = [-1 1];                              
      X_new_1 = zeros(1,D);
      X_new_2 = zeros(1,D);
      
      for d = 1:D
          mut_d = D*(mutant-1)+d;              
          sign1 = A(randi([1,2],1,1));
          X_new_1(d) = X_end(mut_d) + sign1*delta_x(d);
          X_new_2(d) = X_end(mut_d) - sign1*delta_x(d);
      end
      mut_d1 = D*(mutant-1)+1;
      mut_dD = D*mutant;
      X_end(mut_d1 : mut_dD) = X_new_1;
      X_end = [X_end X_new_2];
      N_new = N_end(mutant)/2;
      N_end(mutant) = N_new;
      N_end = [N_end N_new];
      m = m+1;
      mD = m*D;
    end
    
    p = p  + 1
end
 
%% BUILDING A VIDEO
video = VideoWriter ('Equillibration_in_D2.avi', 'Uncompressed AVI'); 
video.FrameRate = 15;      
open(video)                
writeVideo(video, F);       
close(video)               
 

