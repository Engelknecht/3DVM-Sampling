function results = runLoadStep(COORD, SURF, ELEM, NEUMANN, Q, ...
                       elem_type, poisson, traction_force, ...
                       E_elem, Y_elem, a_elem, ...
                       size_xy, size_z, size_hole, sample_id)

results = struct();  % Initialisiere als leeren Struct

% --- Referenzdaten vorbereiten ---
n_n   = size(COORD,2);
n_e   = size(ELEM,2);
[Xi, WF]     = quadrature_volume(elem_type);
[Xi_s, WF_s] = quadrature_surface(elem_type);
[HatP, DHatP1, DHatP2, DHatP3] = local_basis_volume(elem_type, Xi);
[HatP_s, DHatP1_s, DHatP2_s]   = local_basis_surface(elem_type, Xi_s);
n_q   = numel(WF);
n_int = n_e * n_q;

% --- Materialparameter ---
shear_elem = E_elem ./ (2*(1+poisson));
bulk_elem  = E_elem ./ (3*(1-2*poisson));
shear = repelem(shear_elem, n_q);
bulk  = repelem(bulk_elem,  n_q);
a     = repelem(a_elem,     n_q);
Y     = repelem(Y_elem,     n_q);

fprintf('\n[Sample %d] Elemente: %d | Integrationspunkte: %d\n', sample_id, n_e, n_int);

% --- Elastische Steifigkeitsmatrix ---
tic;
[K_elast,B,WEIGHT,iD,jD,D_elast] = elastic_stiffness_matrix(ELEM,COORD,...
                          shear,bulk,DHatP1,DHatP2,DHatP3,WF);
assembly_elast_time = toc;

% --- Neumann-Kräfte ---
n_e_s = size(NEUMANN,2); n_q_s = length(WF_s); n_int_s = n_e_s * n_q_s;
f_t_int = traction_force' * ones(1,n_int_s);
f_t = vector_traction(NEUMANN,COORD,f_t_int,HatP_s,DHatP1_s,DHatP2_s,WF_s);

% --- Initialisierung FEM ---
zeta = [0:0.1:1, 0.9:-0.1:-1, -0.9:0.1:0];
n_step = length(zeta);
alpha = zeros(1,n_step);
U = zeros(3,n_n); dU = U; U_old = U;
F = zeros(3,n_n);
E = zeros(6,n_int);
Ep_old = zeros(6,n_int);
Hard_old = zeros(6,n_int);
assembly = zeros(20*n_step,2);
assembly_step = 0;

% --- Lastschritte ---
for i = 2:n_step
    fprintf('[Sample %d | Step %d]\n', sample_id, i);
    f = zeta(i) * f_t;
    U_it = U; it = 0;

    while true
        tic;
        E(:) = B * U_it(:);
        [S, DS, IND_p] = constitutive_problem(E, Ep_old, Hard_old, shear, bulk, a, Y);
        vD = repmat(WEIGHT, 36, 1) .* DS;
        D_p = sparse(iD(:), jD(:), vD(:), 6*n_int, 6*n_int);
        K_tangent = K_elast + B' * (D_p - D_elast) * B;
        assembly_time = toc;

        n_plast = length(WEIGHT(IND_p));
        assembly_step = assembly_step + 1;
        assembly(assembly_step,:) = [n_plast, assembly_time];

        F(:) = B' * reshape(repmat(WEIGHT,6,1) .* S, 6*n_int, 1);
        dU(Q) = K_tangent(Q,Q) \ (f(Q) - F(Q));
        U_new = U_it + dU;

        q1 = sqrt(dU(:)' * K_elast * dU(:));
        q2 = sqrt(U_it(:)' * K_elast * U_it(:));
        q3 = sqrt(U_new(:)' * K_elast * U_new(:));
        criterion = q1 / (q2 + q3);

        if criterion < 1e-12
            break
        end

        it = it + 1;
        if it > 50
            warning('Sample %d: Newton did not converge.', sample_id);
            results.failed = true;
            return;
        end

        U_it = U_new;
    end

    U_old = U; U = U_it;
    E(:) = B * U(:);
    [S, DS, IND_p, Ep, Hard] = constitutive_problem(E, Ep_old, Hard_old, shear, bulk, a, Y);
    Ep_old = Ep; Hard_old = Hard;
    alpha(i) = f_t(Q)' * U(Q);
end

% --- Ergebnisgrößen speichern ---
stress_vm = sqrt((S(1,:)-S(2,:)).^2 + (S(2,:)-S(3,:)).^2 + (S(3,:)-S(1,:)).^2 ...
             + 6*(S(4,:).^2 + S(5,:).^2 + S(6,:).^2)) / sqrt(2);
strain_eqv = sqrt(sum(E.^2, 1))';

if  mod(i-1,10)==0           %(mod(i,ceil((n_step)/6))==0)
        % displacement
        %U_total = sqrt(U(1,:).^2 + U(2,:).^2 + U(3,:).^2);
        %draw_quantity(COORD,SURF,10*U,U_total,elem_type,size_xy,size_z,size_hole)
               
        % hardening
        if mod(sample_id, 100) == 0
        Hard_node = transformation(sqrt(sum(Hard.^2)),ELEM,WEIGHT); 
        draw_quantity(COORD,SURF,10*U,zeros(size(Hard_node))+Hard_node,elem_type,size_xy,size_z,size_hole) 
        colorbar off; colorbar('location','south')
        end
        
%         fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print('-painters','-dpdf',strcat('figures/fig_VM_3D_hardening_',num2str(i),'_level_',num2str(level)))

    end

results.failed                = false;
results.alpha                 = alpha;
results.zeta                  = zeta;
results.U                     = U;
results.assembly              = assembly;
results.n_int                 = n_int;
results.n_plast               = assembly(:,1);
results.assembly_step         = assembly_step;
results.assembly_elast_time   = assembly_elast_time;
results.stress_vm             = stress_vm(:);
results.strain                = E;
results.strain_eqv            = strain_eqv;

end
