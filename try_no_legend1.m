global D;
global b;
global sigma_comp;
 
data =  load('/home/ane/MathWork/diploma_files/no_legend/HSS/eq/HSS2time681.txt');
b = [0.5129 0.5054; -0.0446 -0.1449];
m = 10;
N_end = data(6:15);
X_end = data(16:35);
x_cc = [0 0];




%b = [0.5129 0.5054; -0.0446 -0.1449];   %HSS
%m = 1;
%N_end = 0.3;
%X_end = [-0.1070 1.4071];
%x_cc = [0 0];

%b = [-0.4796 0.7780; -0.6584 -0.4313];   %LSS
%m = 2;
%N_end = [0.9179 0.9178];
%X_end = [0.5069 0.2021 -0.5057 -0.1860];
%x_cc = [0 0];


% the number of clusters
m_plot = m;
mD = m*D;
sigma_comp = 0.5.*ones(1,D);                           % competition sigma
%x_cc = 0.*ones(1,D);                                   % x_cc - carrying capacity point with D-coordinates
dist_to_join = 0.001;                                   % critical distance to join clusters
dist_to_sep = 0.001;
dist_to_plot = 0.1;
min_size = 0.000001;                                   % critical population size before disappearence
ax = 2;                                                % scale of axis on the plot
axx1= x_cc(1) - 2;
axx2 = x_cc(1) + 2;                                                % scale of axis on the plot
axy1 = x_cc(2) - 2;
axy2 = x_cc(2) + 2;

m_max = 49;
P = 700;
p = 1;
 
delta_x = 0.015;
%goHSS2 0.005
%goHSS2 0.020
%goHSS2 0.050
rr = rand(1);
        cf1 = sqrt(rr)*sign(log(rand(1)+0.5));
        cf2 = sqrt(1 - rr)*sign(log(rand(1)+0.5));
        delta_x1 = cf1*delta_x;
        delta_x2 = cf2*delta_x;
        % x1^2 + x2^2 = delta_x^2
        % delta_x = 0.01, delta_x1^2 + delta_x2^2 = delta_x^2
        % delta_x1^2 = cf1*delta_x^2



speed_ev = zeros((P-1),2);
%X_end = normrnd(0,1,1,mD);                             % initial phenotypes
%N_end = 0.3.*ones(1,m);                                % initial population density of clusters
X_plot = X_end; 
 
% time cycle
tspan_ecology = [0 500];
tspan_evolution = [0 0.1];
 
 
%% FIRST PLOT FIGURE
for k = 1:length(N_end)
        g = 2*k;
        size = 16*N_end(k)+4;
        %size = (17/6)*(log10(N_end(k))+6)+3;         % logarithmic scale!!!!!
        set(0,'DefaultAxesFontSize',16);
        plot (X_end(g-1), X_end(g), 'or', 'MarkerSize', size, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'w');
        axis([axx1 axx2 axy1 axy2]);
        grid on;
        %grid minor;
        hold on;
end
plot (x_cc(1,1), x_cc(1,2), 'or', 'MarkerSize', 6, 'Marker','o','Color','r', 'LineWidth', 1.5);
 
set(gcf,'color','w');            % defines properties of a figure (find figure properties, a lot)
drawnow;                         % Complete pending drawing events
F(1) = getframe(gcf);            % gcf - current figure, save it in vector of frames
pause(.1)


%% TIME STEPS CYCLE
while p < P
    hold off;
    if p > 0
        x_cc(1) = x_cc(1) + delta_x1;
        x_cc(2) = x_cc(2) + delta_x2;
    end
    axx1= x_cc(1) - 2;
    axx2 = x_cc(1) + 2;                                                % scale of axis on the plot
    axy1 = x_cc(2) - 2;
    axy2 = x_cc(2) + 2;
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
        if N_end(q) < min_size           %from what moment we decide to kill the cluster
            N_end(q) = 0;
            X_end( (D*(q-1)+1):q*D ) = 0;
        end
    end
    N_end(N_end == 0) = [];              % delete zeros from N-end vector
    X_end(X_end == 0) = [];
    m = length(N_end);
    mD = m*D;
    X_end_previous = X_end;               % for evolutionary speed 
    
    %% EVOLUTIONARY STEP
    [t X] = ode45(@Phenotypes_cl, tspan_evolution, X_end, [], D, b, m, mD, sigma_comp, N_end, x_cc);
    X_end = X(end, :);
    X = [];
    
    %% JOINING CLOSE CLUSTERS
    if rem(p,10) == 5     % how often we check them to join
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
        m
        m_plot
    
    %% PLOTTING OF THIS CYCLE'S FIGURE
    for k = 1:(length(N_plot))
        if k ==1
            g = 2*k;
            size = 16*N_plot(k)+4;
            %size = (17/6)*(log10(N_end(k))+6)+3;         % logarithmic scale!!!!!
            plot (X_plot(g-1), X_plot(g), 'or', 'MarkerSize', size, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'w');
            set(0,'DefaultAxesFontSize',16);
            axis([axx1 axx2 axy1 axy2]);
            grid on;
            hold on;
        else
            g = 2*k;
            size = 16*N_plot(k)+4;
            %size = (17/6)*(log10(N_end(k))+6)+3;         % logarithmic scale!!!!!
            plot (X_plot(g-1), X_plot(g), 'or', 'MarkerSize', size, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'w');
            set(0,'DefaultAxesFontSize',16);
            axis([axx1 axx2 axy1 axy2]);
            grid on;
            hold on;
        end
    end
    plot (x_cc(1,1), x_cc(1,2), 'or', 'MarkerSize', 6, 'Marker','o', 'Color', 'r', 'LineWidth', 1.5);
    %h = legend([h1 h3],{'clusters','centre of C.C.'});
    %set(h,'FontSize',14,'FontName','Serif');
    xlabel('Phenotype 1','FontSize', 14, 'FontName', 'Serif', 'FontWeight', 'bold');
    ylabel('Phenotype 2','FontSize', 14, 'FontName', 'Serif', 'FontWeight', 'bold');
    %title('Attainment of equilibrium diversity', 'FontName', 'Serif', 'FontWeight', 'normal');
    
    %axis([-ax ax -ax ax])                % axe X from -5 to 5, axe y from -6 to 6
    %grid on;
    %grid minor;
    set(gcf,'color','w');            % defines properties of a figure (find figure properties, a lot)
    drawnow;                         % Complete pending drawing events
    
    if p > 0
        F(p) = getframe(gcf);            % gcf - current figure, save it in vector of frames
        pause(.1)
    end
    %N_end = [0,1] -> b=3, a=17  
    
    %% MUTANT APPEARENCE
    if m < 7^D %& rem(p,10) == 2
      mutant = randi ([1, m], 1, 1);           % choose the cluster to devide
      coef = rand(1,D);                        % coefficients for each dimention
      coef = coef./(sum(coef)/D);              % make some of coefficients equal to 2
      eq_delta_x = dist_to_sep/(2*sqrt(D));
      delta_x = eq_delta_x.*sqrt(coef);   
      A = [-1 1];                              % for direction of coordiate's change
      X_new_1 = zeros(1,D);
      X_new_2 = zeros(1,D);
      
      for d = 1:D
          mut_d = D*(mutant-1)+d;              % coordinates in X_end of cluster ready to mutate
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
    
    tspan_ecology = [0 500];
           
    luc = 1 + D*D + m_plot + D*m_plot;
    output(1, 1:D*D) = (b(:).'); 
    output(1, (D*D+1)) = m_plot;
    output(1, (D*D + 2) : (m_plot + 1 + D*D)) = N_plot;
    output(1, (2 + m_plot + D*D) : (1 + m_plot + D*D + m_plot*D)) = X_plot;
    output = [output x_cc];
    
    if rem(p,20) == 1
        name_file = sprintf('/home/ane/MathWork/diploma_files/no_legend/HSS/0015/go3LSS2time%d.txt',p);
        parsave(name_file,output);
    end
      
    p = p  + 1
end
 
%% BUILDING A VIDEO
video = VideoWriter ('/home/ane/MathWork/diploma_files/no_legend/HSS/0015/go3LSS2_no_legend.avi', 'Uncompressed AVI'); % create a videofile
video.FrameRate = 15;      % the speed of frame's changing
open(video)                % open videofile 
writeVideo(video, F);      % write the video by getting together frames fom F 
close(video)               % close the video 
 
%plot(speed_ev(:,1), speed_ev(:,2), '.')
%means that it was originally under the %
