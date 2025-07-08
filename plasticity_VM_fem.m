%
% The main input data 
%
clc;clear all;

%test
%MI = generate_multiindex(4, 4)

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
M = 5;        % Eigenwerte die betrachtet werden 
lc = 1.0;     % Korrelationslänge
scale_E = 0.05;
scale_Y = 0.05;
scale_a = 0.05;
rho_EY = 0.6; 
rho_Ea = 0.4; 
rho_Ya = 0.5;

%% === PCE Parameter ===
p_order = 2;  %lineares System
Xi_dim = 3*M; %lineares System 
xi_cdf  = @(x) normcdf(x,0,1);
%multi_idx = generate_multiindex(Xi_dim, p_order);
%P = size(multi_idx,1);
%N_samples = 3*P;    

%QoI_samples = zeros(N_samples,1);
%Xi_samples  = zeros(Xi_dim,N_samples);
%
%% Mesh generation

  % the mesh generation depending prescribed finite elements        
  switch(elem_type)
    case 'P1'
        [COORD,ELEM,SURF,NEUMANN,Q]=mesh_P1(level,size_xy,size_z,size_hole);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD,ELEM,SURF,NEUMANN,Q]=mesh_P2(level,size_xy,size_z,size_hole);
        fprintf('P2 elements: \n')
    case 'Q1'
        [COORD,ELEM,SURF,NEUMANN,Q]=mesh_Q1(level,size_xy,size_z,size_hole);
        fprintf('Q1 elements: \n')
    case 'Q2'
        [COORD,ELEM,SURF,NEUMANN,Q]=mesh_Q2(level,size_xy,size_z,size_hole);
        fprintf('Q2 elements: \n')
    otherwise
          disp('bad choice of element type');
  end  

  % Knotenzahl
  n_n = size(COORD,2); % number of nodes
  x = COORD(1,:)';     % x Koordinate 
  y = COORD(2,:)';     % y Koordinate
  points = [x y];

%% KLE

cov_func = @(p1, p2) exp(-norm(p1 - p2)^2 / lc^2); % Gauß
CovMat = zeros(n_n);
for i = 1:n_n
  for j = i:n_n
    val = cov_func(points(i,:), points(j,:));
    CovMat(i,j) = val;
    CovMat(j,i) = val;
  end
end


[Phi, Lambda] = eig(CovMat);
[lambda_sorted, idx] = sort(diag(Lambda), 'descend');
Phi = Phi(:, idx(1:M));
lambda_sorted = lambda_sorted(1:M);

% Beziehung 
Sigma_kle = [eye(M), rho_EY*eye(M), rho_Ea*eye(M);
             rho_EY*eye(M), eye(M), rho_Ya*eye(M);
             rho_Ea*eye(M), rho_Ya*eye(M), eye(M)];

L_kle = chol(Sigma_kle, 'lower');
Z = randn(3*M, 1);
xi_all = L_kle * Z;
xi_E = xi_all(1:M);
xi_Y = xi_all(M+1:2*M);
xi_a = xi_all(2*M+1:end);

% KLE-Felder
E_nodes = zeros(n_n,1);
Y_nodes = zeros(n_n,1);
a_nodes = zeros(n_n,1);

for i = 1:M
  E_nodes = E_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_E(i);
  Y_nodes = Y_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_Y(i);
  a_nodes = a_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_a(i);
end

% Streuung kontrollieren:
% Skalieren auf gewünschte Streuung
%neu

sigma_E = scale_E * mu_E;      % z. B. 0.05·μ
sigma_Y = scale_Y * mu_Y;
sigma_a = scale_a * mu_a;

% Faktor, damit Var(E)=sigma_E² gilt    (λ_sum = Summe der λᵢ)
lambda_sum = sum(lambda_sorted);        % gleiche Summe für alle Knoten,
                                        % weil φs orthonormal bzgl. diskreter Punkte
fac_E = sigma_E / sqrt(lambda_sum);
fac_Y = sigma_Y / sqrt(lambda_sum);
fac_a = sigma_a / sqrt(lambda_sum);

E_nodes = mu_E + fac_E * (Phi * (sqrt(lambda_sorted) .* xi_E));
Y_nodes = mu_Y + fac_Y * (Phi * (sqrt(lambda_sorted) .* xi_Y));
a_nodes = mu_a + fac_a * (Phi * (sqrt(lambda_sorted) .* xi_a));

% bis hier neu

% Auf Elemente interpolieren
n_e = size(ELEM,2);
E_elem = zeros(1, n_e);
Y_elem = zeros(1, n_e);
a_elem = zeros(1, n_e);
for e = 1:n_e
  node_ids = ELEM(:,e);
  E_elem(e) = mean(E_nodes(node_ids));
  Y_elem(e) = mean(Y_nodes(node_ids));
  a_elem(e) = mean(a_nodes(node_ids));
end

%% Postprocessing  KLE Sortieren nach x für Überischt

% Eigenwerte x-Achse
y_vals = unique(round(COORD(2,:), 2));
counts = arrayfun(@(y) sum(abs(COORD(2,:) - y) < 1e-3), y_vals);
[~, max_idx] = max(counts);
line_y = y_vals(max_idx);
line_nodes = abs(COORD(2,:) - line_y) < 1e-3;
x_line = COORD(1, line_nodes);
Phi_line = Phi(line_nodes, 1:M);

% Sortieren nach x für Überischt 
[x_sorted, idx] = sort(x_line);
Phi_sorted = Phi_line(idx, :);
styles = cell(1, M);  % Pre-allokieren
base_styles = {'-','--','-.',':'};
markers = {'o','s','x','^','v','+','*','d','.'};

for i = 1:M
    ls = base_styles{mod(i-1, length(base_styles)) + 1};
    mk = markers{mod(i-1, length(markers)) + 1};
    styles{i} = [ls mk];  % z. B. '--o'
end
figure; hold on;
for i = 1:M
    plot(x_sorted, Phi_sorted(:,i), styles{i}, 'LineWidth', 2, ...
         'DisplayName', sprintf('\\lambda_{%d} = %.3f', i, lambda_sorted(i)));
end
xlabel('x'); ylabel('\phi_i(x)');
legend('show'); title(sprintf('KLE-Eigenfunktionen auf y = %.2f', line_y));
grid on;
%
% Alle Eigenwerte
% Abbildung Eigenwerte KLE Gesamt
figure; hold on;
colors = lines(5); % 5 verschiedene Farbenxx
styles = {'-o','-s','-','--','-.'}; % verschiedene Linienstile

for i = 1:5
    plot(x, Phi(:,i), styles{i}, 'LineWidth', 2, 'DisplayName', ...
        sprintf('\\lambda_{%d} = %.3f', i, lambda_sorted(i)));
end

xlabel('x'); ylabel('\phi_i(x)');
legend('show', 'Location', 'best');
title('Erste 5 Eigenfunktionen der KLE mit zugehörigen Eigenwerten');
grid on;
%
% Visualisierung
figure; scatter3(COORD(1,:), COORD(2,:), Y_nodes, 15, Y_nodes, 'filled');
title('Fließbedingung über Knoten'); colorbar;
figure; scatter3(COORD(1,:), COORD(2,:), a_nodes, 15, a_nodes, 'filled');
title('kinematische Verhärtung über Knoten'); colorbar;
figure; scatter3(COORD(1,:), COORD(2,:), E_nodes, 15, E_nodes, 'filled');
title('E-Modul über Knoten'); colorbar;
% constant tranction on the back side of the body in each direction

%% Data from the reference element
%
  
  % quadrature points and weights for volume and surface integration
  [Xi, WF]     = quadrature_volume(elem_type);
  [Xi_s, WF_s] = quadrature_surface(elem_type);
  
  % local basis functions and their derivatives for volume and surface
  [HatP,DHatP1,DHatP2,DHatP3] = local_basis_volume(elem_type, Xi);
  [HatP_s,DHatP1_s,DHatP2_s]  = local_basis_surface(elem_type, Xi_s);

%
% Number of nodes, elements and integration points + print
%

  n_unknown=length(COORD(Q)); % number of unknowns
  n_e=size(ELEM,2);           % number of elements
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 

  % values of elastic material parameters at integration points 
  shear_elem = E_elem ./ (2*(1+poisson));
  shear = repelem(shear_elem, n_q);
  bulk_elem = E_elem ./ (3*(1-2*poisson));
  bulk = repelem(bulk_elem, n_q);
  a = repelem(a_elem, n_q);
  Y = repelem(Y_elem, n_q);


  fprintf('number of nodes =%d ',n_n);  
  fprintf('\n');   
  fprintf('number of unknowns =%d ',n_unknown);
  fprintf('\n');   
  fprintf('number of elements =%d ',n_e);
  fprintf('\n');   
  fprintf('number of integration points =%d ',n_int);
  fprintf('\n');   

%
% Assembling of the elastic stiffness matrix
%
  
  % values of elastic material parameters at integration points 
  %shear =shear*ones(1,n_int);
  %bulk=bulk*ones(1,n_int);
   
  % stiffness matrix assembly and the assembly time
  tic;     
  [K_elast,B,WEIGHT,iD,jD,D_elast]=elastic_stiffness_matrix(ELEM,COORD,...
                            shear,bulk,DHatP1,DHatP2,DHatP3,WF);  
  assembly_elast_time=toc; 
  fprintf('step =%d ',1);
  fprintf('\n');   
  fprintf('  time spent on K_elast:  %6.1e seconds, ',assembly_elast_time);
  fprintf('\n');   

%
% Assembling of the vector of traction forces
%  
  
  % values of the density function f_t at surface integration points
  n_e_s=size(NEUMANN,2); % number of surface elements
  n_q_s=length(WF_s);    % number of quadrature points on a surface element
  n_int_s=n_e_s*n_q_s ;  % number of integration points on the surface
                         % (on the upper side of the body)
  f_t_int=traction_force'*ones(1,n_int_s); % size(f_V_int)=(3,n_int_s)
  
  % assembling of the vector of traction (surface) forces
  f_t=vector_traction(NEUMANN,COORD,f_t_int,HatP_s,DHatP1_s,DHatP2_s,WF_s);
  
%
% Loading process and Newton's solver
%

  % number of load steps and values of load factors
  d_zeta=1/10;              % constant load increment
  zeta=[0:d_zeta:1, (1-d_zeta):-d_zeta:(-1), (-1+d_zeta):d_zeta:0];
                            % sequence of load factors
  n_step = length(zeta);    % number of load steps
  alpha=zeros(1,n_step);    % work of external forces   
  
  % values of plastic material parematers at integration points
  %a=a*ones(1,n_int);
  %Y=Y*ones(1,n_int);
  
  % initialization
  U = zeros(3,n_n);
  dU = zeros(3,n_n) ; % Newton's increment (in displacement)
  U_old = zeros(3,n_n) ; 
  F = zeros(3,n_n) ;  % vector of internal forces
  E = zeros(6,n_int); % strain tensors at integration points
  Ep_old = zeros(6,n_int); % plastic strain tensors at integration points
  Hard_old = zeros(6,n_int); % hardening tensors at integration points
  
  % storage of assembly time in dependence on plastic integration points
  assembly=zeros(20*n_step,2);
  assembly_step=0;   
  AUX=reshape(1:6*n_int,6,n_int);
  
  % for loop through load steps 
  for i=2:n_step            
      
    fprintf('step =%d ',i);
    fprintf('\n');     
     
    f=zeta(i)*f_t;     % the load vector at step i
    
    % initial displacements
    U_it=U;
      
    % Newton's solver (the semismooth Newton method)
    it=0;              % iteration number
    while true       
        
      % K_consistent tangent stiffness matrix
      tic
      % strain at integration points
      E(:) = B*U_it(:) ;
      % solution of the constitutive problem
      [S,DS,IND_p]=constitutive_problem(E,Ep_old,Hard_old,shear,bulk,a,Y);
      vD = repmat(WEIGHT,36,1).*DS ; 
      D_p = sparse( iD(:),jD(:),vD(:), 6*n_int,6*n_int ) ;   
      K_tangent = K_elast+B'*(D_p-D_elast)*B;          
      assembly_time=toc;
      
      % measuring assembly dependance on plastic integration points
      n_plast=length(WEIGHT(IND_p));
      assembly_step=assembly_step+1;
      assembly(assembly_step,:)=[n_plast assembly_time];      
      
      fprintf('  time spent on K_tangent: %6.1e seconds, ',assembly_time);
      fprintf('  plastic integration points: %d (of %d), ',n_plast,n_int); 
 
      % vector of internal forces
      F(:) = B'*reshape(repmat(WEIGHT,6,1).*S, 6*n_int,1) ;
      
      % Newton's increment
      dU(Q) = K_tangent(Q,Q)\(f(Q)-F(Q)); 
             
      % next iteration
      U_new= U_it + dU ;

      % stopping criterion 
      q1 = sqrt( dU(:)'*K_elast*dU(:) ) ;
      q2 = sqrt(  U_it(:)'*K_elast*U_it(:)  ) ;
      q3 = sqrt( U_new(:)'*K_elast*U_new(:) ) ;
      criterion = q1/(q2+q3);
      
      fprintf('  stopping criterion=%6.1e  ',criterion); 
      fprintf('\n');  
      
      % update of unknown arrays
      U_it=U_new;                                                     
            
      % test on the stopping criterion
      if  criterion < 1e-12
          break
      end
      
      % test on number of iteration
      it=it+1; 
      if  it > 50
          error('The Newton solver does not converge.')
      end         
    end%  true
     
    U_old=U;
    U=U_it;
    E(:) = B*U(:) ;
    [S,DS,IND_p,Ep,Hard]=constitutive_problem(E,Ep_old,Hard_old,shear,bulk,a,Y); 
    Ep_old=Ep;
    Hard_old=Hard;   
    alpha(i)=f_t(Q)'*U(Q);
    
    if  mod(i-1,10)==0           %(mod(i,ceil((n_step)/6))==0)
        % displacement
        %U_total = sqrt(U(1,:).^2 + U(2,:).^2 + U(3,:).^2);
        %draw_quantity(COORD,SURF,10*U,U_total,elem_type,size_xy,size_z,size_hole)
               
        % hardening
        Hard_node = transformation(sqrt(sum(Hard.^2)),ELEM,WEIGHT); 
        draw_quantity(COORD,SURF,10*U,zeros(size(Hard_node))+Hard_node,elem_type,size_xy,size_z,size_hole) 
        colorbar off; colorbar('location','south')
        
%         fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print('-painters','-dpdf',strcat('figures/fig_VM_3D_hardening_',num2str(i),'_level_',num2str(level)))

    end
         
  end %for
  

%  
% Postprocessing - visualization of selected results
%
      
  % mesh visualization
  draw_mesh(COORD,SURF,elem_type)     
  
  % load path
  figure
  plot(alpha,zeta,'x-');
  hold on; plot(alpha([11 21 31 41]),zeta([11 21 31 41]),'ro'); hold off
  axis tight; enlarge_axis(0.1,0.1);
  xlabel('work of loading'); ylabel('loading scale');
  legend('all time steps','visualized time steps','location','northwest')
  

%   fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
%   print('-painters','-dpdf',strcat('figures/fig_VM_3D_hysteresis_level_',num2str(level)));
      

  
  % total displacements + deformed shape
  U_total = sqrt(U(1,:).^2 + U(2,:).^2 + U(3,:).^2);
  draw_quantity(COORD,SURF,10*U,U_total,elem_type,size_xy,size_z,size_hole)
  
  
  
  % dependance of the assembly time on plastic integration points
  figure
  x=assembly(1:assembly_step,1); y=assembly(1:assembly_step,2); 
  x_ext=n_int; y_ext_meas=assembly_elast_time;
  x_r=x(x>0); y_r=y(x>0);  
 
  X = [ones(length(x_r),1) x_r];
  b = X\y_r;  % linear regression, y ~ b(0)+b(1)*x
  x_e=[x; x_ext];      
  y_e=[ones(length(x_e),1) x_e]*b;
  y_ext=[1 x_ext]*b;
  plot(x,y,'bx',x_e,y_e,'b-',x_ext,y_ext,'bo',x_ext,y_ext_meas,'ro');
  legend('plasticity - measurements','plasticity - linear fit',...
  'plasticity - extrapolated value','elasticity','location','northwest')    
  xlabel('number of plastic/all integration points'); ylabel('assembly time');
  axis tight; enlarge_axis(0.1,0.1); 
  
%   fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
%   print('-painters','-dpdf',strcat('figures/fig_VM_3D_assembly_times_',elem_type,'_level_',num2str(level))); 

  
  
  
  


