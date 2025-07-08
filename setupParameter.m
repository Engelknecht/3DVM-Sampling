function [elem_type, level,mu_E, mu_a, mu_Y, poisson, traction_force, size_xy, size_z, size_hole, M, lc, ...
scale_E, scale_Y, scale_a, rho_EY, rho_Ea, rho_Ya, p_order, Xi_dim, xi_cdf  ] = setupParameter()
%% === Hauptparameter ===

  elem_type='Q2'; % the type of finite elements; available choices:              % 'P1', 'P2', 'Q1', 'Q2'
  level=0;        % a nonnegative integer defining mesh density
  
  % values of elastic material parameters
  mu_E = 206900;                         % Young's modulus
  poisson =  0.29;                       % Poisson's ratio
  
  % plastic material parematers
  mu_a = 20000;
  mu_Y = 235 * sqrt(2/3);
  
  % constant tranction on the back side of the body in each direction
  traction_force = [0,200,0];

  % geometrical parameters (choose only integers, size_hole < size_xy)
  size_xy = 10;      % size of the body in direction x and y 
  size_z = 1;        % size of the body in z-direction  
  size_hole = 5;     % size of the hole in the body
                     % the projection of the hole to the xy-plane is square

%% === KLE mit Korrelation zwischen E, Y, a ===
M = 3;        % Eigenwerte die betrachtet werden 
lc = 2.0;     % Korrelationslänge
scale_E = 0.50;
scale_Y = 0.50;
scale_a = 0.50;
rho_EY = 0.6; 
rho_Ea = 0.4; 
rho_Ya = 0.5;

%% === PCE Parameter ===
p_order = 1;  %lineares System
Xi_dim = 3*M; %lineares System 
xi_cdf  = @(x) normcdf(x,0,1);
%multi_idx = generate_multiindex(Xi_dim, p_order);
%P = size(multi_idx,1);
%N_samples = 3*P;    

%QoI_samples = zeros(N_samples,1);
%Xi_samples  = zeros(Xi_dim,N_samples);
%
end

