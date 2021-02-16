% Amplitude list is set manually 
% Constat_adiabaticity_profile.m to calculate it
% or just create Ampl variable and set the proper data
clearvars -except Ampl

%{
  This program is to calculate how resulting magnetization of X nuclei
  depends on some parameters:
    \tau_{sweep}
%}

%% Setting all initial parameters
% Setting spin system
N = 3; I = 1/2;
[Ix, Iy, Iz, ~, ~] = Operators(N,I);

% Setting spin system parameters for fumarate
J_12=15.7; J_13=6.6; J_23=3.2;

% Setting J-coupling Hamiltonian
H0 = J_12*(Ix{1}*Ix{2} + Iy{1}*Iy{2} + Iz{1}*Iz{2}) +...
     J_13*(Ix{1}*Ix{3} + Iy{1}*Iy{3} + Iz{1}*Iz{3}) +...
     J_23*(Ix{2}*Ix{3} + Iy{2}*Iy{3} + Iz{2}*Iz{3});

H1=(1)*(-42.57747892*(Iz{1} + Iz{2}) -10.7084*Iz{3});

% Setting initial density matrix 
ps=[0 0 0 0; 0 1/2 -1/2 0; 0 -1/2 1/2 0; 0 0 0 0];
pcarb=[1/2 0; 0 1/2];
rho_ini=kron(ps,pcarb);

% Diagonalization of the density matrix in the bubbling field
B_0 = Ampl(end);
[eigen_vector_ini, eigen_value_ini] = eig(H0);
rho_ini=inv(eigen_vector_ini)*rho_ini*eigen_vector_ini;
rho_ini=diag(diag(rho_ini));
rho_ini=eigen_vector_ini*rho_ini*inv(eigen_vector_ini);

% Setting the number of point to vary sweep
N_sweep = 100;
% Setting the max value
t_sweep_max = 3;

t_sweep_list = linspace(0, t_sweep_max, N_sweep);

% Setting variable field
I_z = linspace(0, 0, N_sweep)';

for p = 1:N_sweep
    %Initial density matrices
    rho = rho_ini;
    
    % Setting sweeping time
    t_sweep = t_sweep_list(p);
    Nsteps = numel(Ampl);
    dt = t_sweep/(Nsteps-1);
    for q = 1:Nsteps
    H = H0 + H1*Ampl(q);
    rho = expm(-1i*2*pi*dt*H)*rho*expm(1i*2*pi*dt*H);
    end
    
I_z(p) = real(trace(rho*Iz{3}));
end

figure(2)
% The absolute polarization value is negative
% but here I change the sign
plot(t_sweep_list, -I_z);
title('Offset-field dependency')
xlabel('\tau_{sweep} / sec');
ylabel('^{13}Ñ polarization');
