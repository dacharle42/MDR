function [t, x] = mNF_model(t_end_mNF,dt,P_mNF,N_mNF,G_mNF)
%Solves a system of coupled ODEs via MATLAB solver (default: ode45)
opts = odeset('RelTol',1.e-6);
[t, x] = ode15s(@equations,0:dt:t_end_mNF,[P_mNF N_mNF G_mNF],opts);

end

function dx = equations(t,x)
    
    global r_P_N_mNF r_N_P_mNF r_G_N_mNF k_N_mNF g_N_mNF r_G_P_mNF k_G_mNF g_G_mNF N_tot_max_mNF
    
    dx = zeros(3,1);
    
    if (x(1) + x(2) + x(3) < N_tot_max_mNF)
    
        dx(1) = r_P_N_mNF*x(2) - r_N_P_mNF*x(1) - r_G_P_mNF*x(1); %persister cells
        dx(2) = - r_P_N_mNF*x(2) + r_N_P_mNF*x(1) - r_G_N_mNF*x(2) + k_N_mNF*x(2) - g_N_mNF*x(2); %nongenetically drug resistant cells
        dx(3) = r_G_P_mNF*x(1) + r_G_N_mNF*x(2) + k_G_mNF*x(3) - g_G_mNF*x(3); %genetically drug resistant cells
    
    elseif (x(1) + x(2) + x(3) >= N_tot_max_mNF)
        
        dx(1) = 0;
        dx(2) = 0;
        dx(3) = 0;
        
    end
    
end