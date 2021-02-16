% Amplitude list is set manually 
% Constat_adiabaticity_profile.m to calculate it
% or just create Ampl variable and set the proper data
clearvars -except Ampl

%{
  This program is to calculate how resulting magnetization of X nuclei
  depends on some parameters:
    z-field offset
    transverse field
  To choose what to vary, set case_var respetively:
    B_z
    B_x
%}
case_var = 'B_z';

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

%%Setting case and calculating outcome
switch case_var
    case 'B_z'
        % Zeeman variable Hamiltonain 
        H_var = (1)*(-42.57747892*(Iz{1}+Iz{2})-10.7084*Iz{3});
    case 'B_x'
        % Zeeman variable Hamiltonain 
        H_var = (1)*(-42.57747892*(Ix{1}+Ix{2})-10.7084*Ix{3});
end
% Setting the number of point to vart field
N_def = 100;
        
% Setting sweeping time
t_sweep = 3;
Nsteps = numel(Ampl);
dt = t_sweep/(Nsteps-1);

% Setting variable field
B = linspace(-0.2,0.2,N_def)';
I_z = linspace(0,0,N_def)';

for p = 1:N_def
    %Initial density matrices
    rho = rho_ini;    
    for q = 1:Nsteps
    H = H0 + H1*Ampl(q) + H_var*B(p);
    rho = expm(-1i*2*pi*dt*H)*rho*expm(1i*2*pi*dt*H);
    end
    
I_z(p) = real(trace(rho*Iz{3}));
end

figure(2)
% The absolute polarization value is negative
% but here I change the sign
plot(B, -I_z);
title('Offset-field dependency')
xlabel('field offset / \muT');
ylabel('^{13}Ñ polarization');
