function [t, x] = mPF_model(t_end_mPF,dt,P_mPF,N_mPF,G_mPF)
%Solves a system of coupled ODEs via MATLAB solver (default: ode45)
opts = odeset('RelTol',1.e-6);
[t, x] = ode15s(@equations,0:dt:t_end_mPF,[P_mPF N_mPF G_mPF],opts);

end

function dx = equations(t,x)
    
    global r_P_N_mPF r_N_P_mPF r_G_N_mPF k_N_mPF g_N_mPF r_G_P_mPF k_G_mPF g_G_mPF N_tot_max_mPF
    
    dx = zeros(3,1);
    
    if (x(1) + x(2) + x(3) < N_tot_max_mPF)
     

        dx(1) = r_P_N_mPF*x(2) - r_N_P_mPF*x(1) - r_G_P_mPF*x(1); %persister cells
        dx(2) = - r_P_N_mPF*x(2) + r_N_P_mPF*x(1) - r_G_N_mPF*x(2) + k_N_mPF*x(2) - g_N_mPF*x(2); %nongenetically drug resistant cells
        dx(3) = r_G_P_mPF*x(1) + r_G_N_mPF*x(2) + k_G_mPF*x(3) - g_G_mPF*x(3); %genetically drug resistant cells
        
    elseif (x(1) + x(2) + x(3) >= N_tot_max_mPF)
        
         dx(1) = 0;
         dx(2) = 0;
         dx(3) = 0;
        
    end
    
end