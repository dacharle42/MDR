%% Daniel Charlebois - Stochastic Population Dynamics Model - Matlab vR2018b 
%Stochastic population algorithm - numerically simulates the ordinary differential equation model describing the growth and evolution 
%of the mammalian Negative Feedback (mNF) and mammalian Positive Feedback (mPF) gene circuits in Chinease Hamster Ovary (CHO) cells 
%under various Puromycin conditions. This script calls the functions "mPF_model.m" and "mNF_model.m".

clc; clear all; close all;

tic

%% parameters
%global
global r_P_N_mPF r_N_P_mPF r_G_N_mPF k_N_mPF g_N_mPF r_G_P_mPF k_G_mPF g_G_mPF N_tot_max_mPF %mPF
global r_P_N_mNF r_N_P_mNF r_G_N_mNF k_N_mNF g_N_mNF r_G_P_mNF k_G_mNF g_G_mNF N_tot_max_mNF %mNF
%local
Puro = 0; %[0 10 22.5 35 50]ug/ml
td = 18; %cell division time (hr)
mu = 1/td; %growth rate (1/hr)
dt = 1; %(hr)
num_runs = 6;
N_tot_max_data_mPF(num_runs) = 0;
N_tot_max_data_mNF(num_runs) = 0;

for k = 1:num_runs
    %parameter values and initial conditions for each Puromycin condition
    if Puro == 0
        %noise terms
        sigma = 100; sigma2 = 100;
        
        %mPF
        r_P_N_mPF = 10^(-4); r_N_P_mPF = 10^(-2); r_G_N_mPF = 0; k_N_mPF = mu; g_N_mPF = 1*10^(-2); r_G_P_mPF = 0; k_G_mPF = mu; g_G_mPF = 1*10^(-2);     
        %carrying capacity
        N_tot_max_mPF = 6*10^4+randn*sigma2; 
        %experiment duration (hr)
        t_end_mPF = 100; t_data_mPF=0:dt:t_end_mPF;
        %initial fraction and number of cells surviving initial Puromycin treatment
        N_tot_init_mPF = 0.1*10^4; frac_survive_mPF = 1; N_survive_mPF = round((frac_survive_mPF*N_tot_init_mPF)+(randn*sigma)); 
        if (N_survive_mPF < 0), N_survive_mPF = 0; end %ensure a non-negative number of surviving cells at any level of sigma
        %initial D (dead), P (persister), N (nongenetically resistant), and G (genetically resistant) fractions (A)/numbers of cells
        AD_mPF = 1 - frac_survive_mPF;
        AP_mPF = 0; P_mPF = round(AP_mPF*N_survive_mPF); 
        AG_mPF = 0; G_mPF = 0; %clonal cell line with no pre-existing drug resistance mutations
        AN_mPF = 1-AP_mPF-AG_mPF; N_mPF = round(AN_mPF*N_survive_mPF);
        
        %mNF
        r_P_N_mNF = 10^(-4); r_N_P_mNF = 10^(-2); r_G_N_mNF = 0; k_N_mNF = mu; g_N_mNF = 1*10^(-2); r_G_P_mNF = 0; k_G_mNF = mu; g_G_mNF = 1*10^(-2);    
        N_tot_max_mNF = 6*10^4+randn*sigma2;
        t_end_mNF = 100; t_data_mNF=0:dt:t_end_mNF;
        N_tot_init_mNF = 0.085*10^4; frac_survive_mNF = 1; N_survive_mNF = round((frac_survive_mNF*N_tot_init_mNF)+(randn*sigma)); 
        if (N_survive_mNF < 0), N_survive_mNF = 0; end
        AD_mNF = 1 - frac_survive_mNF;
        AP_mNF = 0; P_mNF = round(AP_mNF*N_survive_mNF); 
        AG_mNF = 0; G_mNF = 0;
        AN_mNF = 1-AP_mNF-AG_mNF; N_mNF = round(AN_mNF*N_survive_mNF);
    elseif Puro == 10       
        sigma = 500; sigma2 = 500;
        
        %mPF
        r_P_N_mPF = 10^(-4); r_N_P_mPF = 10^(-2); r_G_N_mPF = 10^(-7); k_N_mPF = 0.825*mu; g_N_mPF = 2.5*10^(-2); r_G_P_mPF = 5*10^(-6); k_G_mPF = 0.825*mu; g_G_mPF = 1*10^(-2); 
        N_tot_max_mPF = 3*10^4+randn*sigma2;
        t_end_mPF = 255; t_data_mPF=0:dt:t_end_mPF; 
        N_tot_init_mPF = 0.15*10^4; frac_survive_mPF = 0.95; N_survive_mPF = round((frac_survive_mPF*N_tot_init_mPF)+(randn*sigma)); 
        if (N_survive_mPF < 0), N_survive_mPF = 0; end
        AD_mPF = 1 - frac_survive_mPF;
        AP_mPF = 0; P_mPF = round(AP_mPF*N_survive_mPF); 
        AG_mPF = 0; G_mPF = 0;
        AN_mPF = 1-AP_mPF-AG_mPF; N_mPF = round(AN_mPF*N_survive_mPF);
                
        %mNF
        r_P_N_mNF = 10^(-4); r_N_P_mNF = 10^(-2); r_G_N_mNF = 10^(-6); k_N_mNF = 0.825*mu; g_N_mNF = 2.5*10^(-2); r_G_P_mNF = 5*10^(-6); k_G_mNF = 0.825*mu; g_G_mNF = 1*10^(-2); 
        N_tot_max_mNF = 4*10^4+randn*sigma2;        
        t_end_mNF = 200; t_data_mNF=0:dt:t_end_mNF;
        N_tot_init_mNF = 0.2*10^4; frac_survive_mNF = 0.99; N_survive_mNF = round((frac_survive_mNF*N_tot_init_mNF)+(randn*sigma));
        if (N_survive_mNF < 0), N_survive_mNF = 0; end
        AD_mNF = 1 - frac_survive_mNF;
        AP_mNF = 0; P_mNF = round(AP_mNF*N_survive_mNF); 
        AG_mNF = 0; G_mNF = 0;
        AN_mNF = 1-AP_mNF-AG_mNF; N_mNF = round(AN_mNF*N_survive_mNF);
        
    elseif Puro == 22.5    
        sigma = 400; sigma2 = 250;
        
        %mPF
        r_P_N_mPF = 10^(-4); r_N_P_mPF = 10^(-2); r_G_N_mPF = 10^(-7); k_N_mPF = 0.6*mu; g_N_mPF = 3*10^(-2); r_G_P_mPF = 5*10^(-6); k_G_mPF = 0.8*mu; g_G_mPF = 1*10^(-2);
        N_tot_max_mPF = 3.25*10^4+randn*sigma2;
        t_end_mPF = 380; t_data_mPF=0:dt:t_end_mPF;
        N_tot_init_mPF = 0.2*10^4; frac_survive_mPF = 0.59; N_survive_mPF = round((frac_survive_mPF*N_tot_init_mPF)+(randn*sigma)); 
        if (N_survive_mPF < 0), N_survive_mPF = 0; end
        AD_mPF = 1 - frac_survive_mPF;
        AP_mPF = frac_survive_mPF; P_mPF = round(AP_mPF*N_survive_mPF); 
        AG_mPF = 0; G_mPF = 0;
        AN_mPF = 1-AP_mPF-AG_mPF; N_mPF = round(AN_mPF*N_survive_mPF);
        
        %mNF
        r_P_N_mNF = 10^(-4); r_N_P_mNF = 10^(-2); r_G_N_mNF = 10^(-6); k_N_mNF = 0.8*mu; g_N_mNF = 3*10^(-2); r_G_P_mNF = 5*10^(-6); k_G_mNF = 0.8*mu; g_G_mNF = 1*10^(-2);
        N_tot_max_mNF = 3*10^4+randn*sigma2;
        t_end_mNF = 225; t_data_mNF=0:dt:t_end_mNF;
        N_tot_init_mNF = 0.25*10^4; frac_survive_mNF = 0.90; N_survive_mNF = round((frac_survive_mNF*N_tot_init_mNF)+(randn*sigma));
        if (N_survive_mNF < 0), N_survive_mNF = 0; end
        AD_mNF = 1 - frac_survive_mNF;
        AP_mNF = 0.8133; P_mNF = round(AP_mNF*N_survive_mNF); 
        AG_mNF = 0; G_mNF = 0;
        AN_mNF = 1-AP_mNF-AG_mNF; N_mNF = round(AN_mNF*N_survive_mNF);
    elseif Puro == 35   
        sigma = 1000; sigma2 = 1000; 
               
        %mPF
        r_P_N_mPF = 10^(-4); r_N_P_mPF = 10^(-2); r_G_N_mPF = 10^(-7); k_N_mPF = 0.45*mu; g_N_mPF = 3.5*10^(-2); r_G_P_mPF = 5*10^(-6); k_G_mPF = 0.45*mu; g_G_mPF = 1*10^(-2);       
        N_tot_max_mPF = 2.5*10^4+randn*sigma2; N_tot_max_data_mPF(k) = N_tot_max_mPF;
        t_end_mPF = 1500; t_data_mPF=0:dt:t_end_mPF;
        N_tot_init_mPF = 0.1*10^4; frac_survive_mPF = 0.37; N_survive_mPF = round((frac_survive_mPF*N_tot_init_mPF)+(randn*sigma)); 
        if (N_survive_mPF < 0), N_survive_mPF = 0; end
        AD_mPF = 1 - frac_survive_mPF;
        AP_mPF = frac_survive_mPF; P_mPF = round(AP_mPF*N_survive_mPF);
        AG_mPF = 0; G_mPF = 0;
        AN_mPF = 1-AP_mPF-AG_mPF ; N_mPF = round(AN_mPF*N_survive_mPF);
        
        %mNF
        r_P_N_mNF = 10^(-4); r_N_P_mNF = 10^(-2); r_G_N_mNF = 10^(-6); k_N_mNF = 0.45*mu; g_N_mNF = 3.5*10^(-2); r_G_P_mNF = 5*10^(-6); k_G_mNF = 0.45*mu; g_G_mNF = 1*10^(-2);      
        N_tot_max_mNF = 2.25*10^4+randn*sigma2; N_tot_max_data_mNF(k) = N_tot_max_mNF;
        t_end_mNF = 1005; t_data_mNF=0:dt:t_end_mNF;
        N_tot_init_mNF = 0.1*10^4; frac_survive_mNF = 0.51; N_survive_mNF = round((frac_survive_mNF*N_tot_init_mNF)+(randn*sigma));
        if (N_survive_mNF < 0), N_survive_mNF = 0; end
        AD_mNF = 1 - frac_survive_mNF;
        AP_mNF = frac_survive_mPF; P_mNF = round(AP_mNF*N_survive_mNF); 
        AG_mNF = 0; G_mNF = 0;
        AN_mNF = 1-AP_mNF-AG_mNF; N_mNF = round(AN_mNF*N_survive_mNF);
    elseif Puro == 50
        sigma = 1000; sigma2 = 500;
        
        %mPF
        r_P_N_mPF = 10^(-4); r_N_P_mPF = 10^(-2); r_G_N_mPF = 10^(-7); k_N_mPF = 0.4*mu; g_N_mPF = 5*10^(-2); r_G_P_mPF = 5*10^(-6); k_G_mPF = 0.4*mu; g_G_mPF = 1*10^(-2); 
        N_tot_max_mPF = 2*10^4+randn*sigma2; N_tot_max_data_mPF(k) = N_tot_max_mPF; 
        t_end_mPF = 1500; t_data_mPF=0:dt:t_end_mPF;
        N_tot_init_mPF = 0.1*10^4; frac_survive_mPF = 0.1; N_survive_mPF = round((frac_survive_mPF*N_tot_init_mPF)+(randn*sigma)); 
        if (N_survive_mPF < 0), N_survive_mPF = 0; end
        AD_mPF = 1 - frac_survive_mPF;
        AP_mPF = frac_survive_mPF; P_mPF = round(AP_mPF*N_survive_mPF);
        AG_mPF = 0; G_mPF = 0;
        AN_mPF = 1-AP_mPF-AG_mPF; N_mPF = round(AN_mPF*N_survive_mPF);
        
        %mNF
        r_P_N_mNF = 10^(-4); r_N_P_mNF = 10^(-2); r_G_N_mNF = 10^(-6); k_N_mNF = 0*mu; g_N_mNF = 5*10^(-2); r_G_P_mNF = 5*10^(-6); k_G_mNF = 0.4*mu; g_G_mNF = 1*10^(-2); 
        N_tot_max_mNF = 2*10^4+randn*sigma2; N_tot_max_data_mNF(k) = N_tot_max_mNF; 
        t_end_mNF = 1500; t_data_mNF=0:dt:t_end_mNF;
        N_tot_init_mNF = 0.1*10^4; frac_survive_mNF = 0; N_survive_mNF = 0;
        if (N_survive_mNF < 0), N_survive_mNF = 0; end
        AD_mNF = 1 - frac_survive_mNF;
        AP_mNF = frac_survive_mPF; P_mNF = round(AP_mNF*N_survive_mNF);
        AG_mNF = 0; G_mNF = 0;
        AN_mNF = 1-AP_mNF-AG_mNF; N_mNF = round(AN_mNF*N_survive_mNF); 
    end

    %% numerical simulation
    %call to ODEsolver - mPF
    [t, x] = mPF_model(t_end_mPF,dt,P_mPF,N_mPF,G_mPF);
    %store data - mPF
    P_data_mPF(:,k) = x(:,1); N_data_mPF(:,k) = x(:,2); G_data_mPF(:,k) = x(:,3);    
    N_tot_data_mPF(:,k) = x(:,1) + x(:,2) + x(:,3);
    %call to ODEsolver - mNF
    [t, x] = mNF_model(t_end_mNF,dt,P_mNF,N_mNF,G_mNF);
    %store data - mNF
    P_data_mNF(:,k) = x(:,1); N_data_mNF(:,k)= x(:,2); G_data_mNF(:,k)= x(:,3);
    N_tot_data_mNF(:,k) = x(:,1) + x(:,2) + x(:,3);
    
end

%% process data 
%Note CHO experiments end when replicate reaches confluency in Puromycin 35
%ug/ml and 50 ug/ml conditions - for all other Puromycin conditions CHO experiments 
%for all replicates were ended at the same timepoint (no processing required).
if Puro == 35 || Puro == 50
    for i = 1:num_runs
        %save array and time value only up to N_max
        %mPF
        idx_hld_mPF = find(N_tot_data_mPF(:,i) >= N_tot_max_data_mPF(i));
        if isempty(idx_hld_mPF) == 1
            idx_mPF(i) = t_end_mPF+1;
        else
            idx_mPF(i) = idx_hld_mPF(1);
        end
        %mNF
        idx_hld_mNF = find(N_tot_data_mNF(:,i) >= N_tot_max_data_mNF(i));
        if isempty(idx_hld_mNF) == 1
            idx_mNF(i) = t_end_mNF+1;
        else
            idx_mNF(i) = idx_hld_mNF(1);
        end
    end
end
    
%% save data
%mPF
save('mPF_pop_model_Puro0.mat','t_data_mPF','P_data_mPF','N_data_mPF','G_data_mPF','N_tot_data_mPF','Puro');
%mNF
save('mNF_pop_model_Puro0.mat','t_data_mNF','P_data_mNF','N_data_mNF','G_data_mNF','N_tot_data_mNF','Puro');

%% plot results
subplot(1,3,1)
if Puro == 35 || Puro == 50
    hold on
    plot(t_data_mPF(1:idx_mPF(1)),P_data_mPF(1:idx_mPF(1),1),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
    plot(t_data_mPF(1:idx_mPF(1)),N_data_mPF(1:idx_mPF(1),1),':','Color',[0 0.7 0],'LineWidth',3) 
    plot(t_data_mPF(1:idx_mPF(1)),G_data_mPF(1:idx_mPF(1),1),'k:','LineWidth',3) 

    for i=2:num_runs
        plot(t_data_mPF(1:idx_mPF(i)),P_data_mPF(1:idx_mPF(i),i),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
        plot(t_data_mPF(1:idx_mPF(i)),N_data_mPF(1:idx_mPF(i),i),':','Color',[0 0.7 0],'LineWidth',3) 
        plot(t_data_mPF(1:idx_mPF(i)),G_data_mPF(1:idx_mPF(i),i),'k:','LineWidth',3) 
    end
    hold off
elseif Puro == 0 || Puro == 10 || Puro == 22.5
    hold on
    plot(t_data_mPF,P_data_mPF(:,1),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
    plot(t_data_mPF,N_data_mPF(:,1),':','Color',[0 0.7 0],'LineWidth',3) 
    plot(t_data_mPF,G_data_mPF(:,1),'k:','LineWidth',3) 

    plot(t_data_mPF,P_data_mPF(:,2:num_runs),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
    plot(t_data_mPF,N_data_mPF(:,2:num_runs),':','Color',[0 0.7 0],'LineWidth',3) 
    plot(t_data_mPF,G_data_mPF(:,2:num_runs),'k:','LineWidth',3) 
    hold off
end
xlabel('time (hours)'); ylabel('cell count');
legend('P', 'N', 'G','Location','northwest')
title('mPF Subpopulation Growth')
subplot(1,3,2)
if Puro == 35 || Puro == 50
    hold on
    plot(t_data_mNF(1:idx_mNF(1)),P_data_mNF(1:idx_mNF(1),1),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
    plot(t_data_mNF(1:idx_mNF(1)),N_data_mNF(1:idx_mNF(1),1),':','Color',[0 0.7 0],'LineWidth',3) 
    plot(t_data_mNF(1:idx_mNF(1)),G_data_mNF(1:idx_mNF(1),1),'k:','LineWidth',3) 

    for i=2:num_runs
        plot(t_data_mNF(1:idx_mNF(i)),P_data_mNF(1:idx_mNF(i),i),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
        plot(t_data_mNF(1:idx_mNF(i)),N_data_mNF(1:idx_mNF(i),i),':','Color',[0 0.7 0],'LineWidth',3) 
        plot(t_data_mNF(1:idx_mNF(i)),G_data_mNF(1:idx_mNF(i),i),'k:','LineWidth',3) 
    end
    hold off
elseif Puro == 0 || Puro == 10 || Puro == 22.5
    hold on
    plot(t_data_mNF,P_data_mNF(:,1),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
    plot(t_data_mNF,N_data_mNF(:,1),':','Color',[0 0.7 0],'LineWidth',3) 
    plot(t_data_mNF,G_data_mNF(:,1),'k:','LineWidth',3) 

    plot(t_data_mNF,P_data_mNF(:,2:num_runs),':','Color',[0.5 0.5 0.5],'LineWidth',3) 
    plot(t_data_mNF,N_data_mNF(:,2:num_runs),':','Color',[0 0.7 0],'LineWidth',3) 
    plot(t_data_mNF,G_data_mNF(:,2:num_runs),'k:','LineWidth',3) 
    hold off
end
xlabel('time (hours)'); ylabel('cell count');
legend('P', 'N', 'G','Location','northwest')
title('mNF Subpopulation Growth')

subplot(1,3,3)
if Puro == 35 || Puro == 50
     hold on
     
       plot(t_data_mPF(1:idx_mPF(1)),N_tot_data_mPF(1:idx_mPF(1),1),'r:','LineWidth',3)
       plot(t_data_mNF(1:idx_mNF(1)),N_tot_data_mNF(1:idx_mNF(1),1),'b:','LineWidth',3)
       
     for i=2:num_runs
       plot(t_data_mPF(1:idx_mPF(i)),N_tot_data_mPF(1:idx_mPF(i),i),'r:','LineWidth',3)
       plot(t_data_mNF(1:idx_mNF(i)),N_tot_data_mNF(1:idx_mNF(i),i),'b:','LineWidth',3)
     end
     hold off
elseif Puro == 0 || Puro == 10 || Puro == 22.5
    hold on
    plot(t_data_mPF,N_tot_data_mPF(:,1),'r:','LineWidth',3)
    plot(t_data_mNF,N_tot_data_mNF(:,1),'b:','LineWidth',3)
 
    plot(t_data_mPF,N_tot_data_mPF(:,2:num_runs),'r:','LineWidth',3)
    plot(t_data_mNF,N_tot_data_mNF(:,2:num_runs),'b:','LineWidth',3)
    hold off
end
xlabel('time (hours)'); ylabel('cell count');
title('mPF & mNF Population Growth')
legend('mPF', 'mNF','Location','northwest')

toc