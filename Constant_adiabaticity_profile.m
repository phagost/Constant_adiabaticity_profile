close all; clear all;

%{
  This program is to calculate constant adiabaticity profile for HP methods
  called field sweep either from or through zero (FS-f-Z or FS-t-Z)
%}

%{
  To chose the method switch the case_var as:
    "FS-f-Z"
    or
    "FS-t-Z"
  And set appropriate starting value ( Ampl(1) )
%}
case_var = 'FS-t-Z'; 

%% Setting all initial parameters
% Setting spin system
N = 3; I = 1/2;
[Ix, Iy, Iz, ~, ~] = Operators(N, I);

% Setting spin system parameters for fumarate
J_12 = 15.7; J_13 = 6.6; J_23 = 3.2;
gyro_H = 42.57747892; gyro_C = 10.7084; % MHz T^{-1}

% Setting J-coupling Hamiltonian
H0 = J_12*(Ix{1}*Ix{2} + Iy{1}*Iy{2} + Iz{1}*Iz{2}) +...
     J_13*(Ix{1}*Ix{3} + Iy{1}*Iy{3} + Iz{1}*Iz{3}) +...
     J_23*(Ix{2}*Ix{3} + Iy{2}*Iy{3} + Iz{2}*Iz{3}); 
 
H1=(1)*(-42.57747892*(Iz{1} + Iz{2}) -10.7084*Iz{3});  
 
% Setting unitary matrix to make Hamiltonian block-diagonal
% i.e. the angular projection is in order {3/2, 1/2, -1/2, -3/2}
U_block = diag(diag(ones(8)));
U_block(4, 4) = 0; U_block(4, 5) = 1;
U_block(5, 5) = 0; U_block(5, 4) = 1;

% change the Hamiltonian basis
H0_new = inv(U_block)*H0*U_block;
H1_new = inv(U_block)*H1*U_block;

% amplitude list in uT
Ampl = zeros(1000,1);     

%{
  The maximal value untill the profile is caluclated
  this should be much higher then the LAC field
  the coordinate of the LAC field may be found here
  https://doi.org/10.1063/1.5089486
%}
B0_max = 2;

%{
  The initial value of the field
  IT DEPENDS ON THE METHODS:
    FS-f-Z:
        zero or near-zero (50 nT in our case)
    FS-t-Z:
        the final field with reverse sign (-B1m)
%}
%Ampl(1) = - B0_max;
Ampl(1) = 0;

% cycle counter
inner_counter = 0;

% time increment interval
dt = 1/20000;


%% The cycle to calculate profile
while true
    inner_counter = inner_counter+1;
    
    % setting the hamiltonian for this iteration
    HH = H0_new + H1_new*Ampl(inner_counter);
    
    % Chosing two central blocks 3x3
    [eigen_vector_1, eigen_value_1] = eig(HH(2:4, 2:4));
    [eigen_vector_2, eigen_value_2] = eig(HH(5:7, 5:7));
    
    % Operator that gives adiabatic paremeter, i.e.
    % dH/dt=-Ixx*dv1(t)/dt
    Op_1 = inv(eigen_vector_1)*(H1(2:4, 2:4))*eigen_vector_1;    
    Op_2 = inv(eigen_vector_2)*(H1(5:7, 5:7))*eigen_vector_2;    
    
    % Calculation of the adiabaticity parameter
    ksi = 0;
    for i = 1:3
        for j = 1:3
            if i == j 
            %if i==j || (i==2&&j==3) || (i==3&&j==2)
                ksi1 = 0;   % this term is null cause Ixx is in its own 
                ksi2 = 0;   % eigenbasis Ixx|i>~|i>
            elseif i ~= j
                % Adiabaticity parameter by def
                switch case_var
                case "FS-t-Z"
                    if i==1 | j==1
                        ksi1 = abs(Op_1(i,j))*abs(Op_1(j,i))/...
                            (eigen_value_1(i,i)-eigen_value_1(j,j))^4;
                    end
                    if i==2 | j==2
                        ksi2 = abs(Op_2(i,j))*abs(Op_2(j,i))/...
                            (eigen_value_2(i,i)-eigen_value_2(j,j))^4;
                    end
                 case "FS-f-Z"
                    if i==1 | j==1
                        ksi1 = abs(Op_1(i,j))*abs(Op_1(j,i))/...
                            (eigen_value_1(i,i)-eigen_value_1(j,j))^4;
                        ksi2 = abs(Op_2(i,j))*abs(Op_2(j,i))/...
                            (eigen_value_2(i,i)-eigen_value_2(j,j))^4;
                    end
                end
            end
            ksi = ksi + ksi1 + ksi2;
        end
    end
    ksi = sqrt(ksi);
    
% Here the number of the amplitude element list is increased
% if the default number is exceeded
    if inner_counter > numel(Ampl)
        u = 2*numel(Ampl);
        Ampldumm = zeros(u,1);
        for l=1:numel(Ampl)
            Ampldumm(l)=Ampl(l);
        end
        Ampl=Ampldumm;
        clearvars Ampldumm
    end

% Calclulating the next B0 value
% or stopping
    if Ampl(inner_counter )< B0_max
        Ampl(inner_counter+1) = Ampl(inner_counter)+1*(1/ksi)*dt;
    elseif Ampl(inner_counter) > B0_max
        Ampl(inner_counter) = B0_max;
        break
    end
end

% Deleting all zero elements
Ampl = Ampl(1:inner_counter);
    
% Interpolation
B1m = Ampl(end);
N_interpol = 2000;
data1 = interp1(linspace(0, 1, numel(Ampl)), Ampl,...
                linspace(0, 1, N_interpol), 'spline');
Ampl = data1';

% Plotting and comparing with linear
Ampllin=linspace(Ampl(1), B1m, N_interpol)';
figure(1);
plot(linspace(0, 1, N_interpol), Ampl,...
     linspace(0, 1, N_interpol), Ampllin)
 title('Constant adiabaticity field profile');
 xlabel('t / \tau_{sweep}');
 ylabel('Magnetic Field / \muT');

%{
Ntimes=600; %The number of checking switching times
%Tvar=linspace(0,TT*15,Ntimes)';
Tvar=linspace(0,6,Ntimes)';
%Arrays for carbon magnetization after sweeping
IzCopt=zeros(Ntimes,1);
IzClin=zeros(Ntimes,1);

%Initial density matrix - singlet between protons and no polarization in
%the 13 carbon
r0=zeros(8);
r0(3,3)=1/4; r0(4,4)=1/4; r0(5,5)=1/4; r0(6,6)=1/4;
r0(3,4)=-1/4; r0(4,3)=-1/4; r0(5,6)=-1/4; r0(6,5)=-1/4;
%
%Change a dens a little bit
[eigen_vector_fin,eigen_value_fin]=eig(H0_new+H1_new*Ampl(end));
r0=inv(eigen_vector_fin)*r0*eigen_vector_fin;
r0=diag(diag(r0));
r0=eigen_vector_fin*r0*inv(eigen_vector_fin);

%Projector
%
Iz_new=inv(U_block)*Iz{3}*U_block;
for q=1:Ntimes
    %Initial density matrices
    rrropt=r0;  rrrlin=r0;
    
    %dummy time for density matrix calculating
%Tdum=linspace(0,Tvar(q),iss-1);
%    dt=handles.Tvar(q)/(inner_counter-1);
    dt=Tvar(q)/(Nsteps-1);
    for is=1:Nsteps-1
    nopt=Ampl(is);
    nlin=Ampllin(is);
    Hopt=H0_new+H1_new*nopt;
    Hlin=H0_new+H1_new*nlin;
    rrropt=expm(-1i*2*pi*dt*Hopt)*rrropt*expm(1i*2*pi*dt*Hopt);
    rrrlin=expm(-1i*2*pi*dt*Hlin)*rrrlin*expm(1i*2*pi*dt*Hlin);
    end
    
%IzCopt(q)=real(trace(rrropt*Iz{3}));
%IzClin(q)=real(trace(rrrlin*Iz{3}));

IzCopt(q)=real(trace(rrropt*Iz_new));
IzClin(q)=real(trace(rrrlin*Iz_new));

end

figure(2)
plot(Tvar,IzCopt,'b',Tvar,IzClin,'r');

%}

%% Here I plot Iz from time dependency
%{
time_switch=3*0.1465;
rrropt=r0;  rrrlin=r0;
    dt=time_switch/(Nsteps-1);
    for is=1:Nsteps
    nopt=Ampl(is)+2;
    nlin=Ampllin(is)+2;
    Hopt=handles.H0+handles.H1*nopt;
    Hlin=handles.H0+handles.H1*nlin;
    rrropt=expm(-1i*2*pi*dt*Hopt)*rrropt*expm(1i*2*pi*dt*Hopt);
    rrrlin=expm(-1i*2*pi*dt*Hlin)*rrrlin*expm(1i*2*pi*dt*Hlin);
    IzCopt_time(is)=real(trace(rrropt*Iz{3}));
    IzClin_time(is)=real(trace(rrrlin*Iz{3}));
    end
    
figure(3)
switching_time=linspace(0,time_switch,Nsteps)';
plot(switching_time,IzCopt_time,'b',switching_time,IzClin_time,'r');
%}