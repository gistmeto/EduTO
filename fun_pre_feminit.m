function [fem] = fun_pre_feminit(inputs,msh)

% Separation of FEM domain and exterior domain
    ele_ind_fem=find(msh.IX(:,5) <=2);
    node_ind_fem= unique([msh.IX(ele_ind_fem,1);msh.IX(ele_ind_fem,2);...
                          msh.IX(ele_ind_fem,3);msh.IX(ele_ind_fem,4)]);
    for i=1:size(node_ind_fem,1); node_ind_inv(node_ind_fem(i))=i; end                     
    fem.extX=msh.X(setdiff([1:size(msh.X,1)]',node_ind_fem),:);

% Setting of FEM variables
    fem.IX=msh.IX(ele_ind_fem,:);       % Element connectivity
    fem.IX=[node_ind_inv(fem.IX(:,[1:4])) fem.IX(:,5)];
    fem.X=msh.X(node_ind_fem,:);        % Nodal coordinates
    fem.nn = size(fem.X,1);                 % Number of node
    fem.ndof = fem.nn*3;                    % Number of degrees of freedom; u, v
    fem.ne = size(fem.IX,1);                % Number of elements    
    fem.edof = [fem.IX(:,1)*3-2 fem.IX(:,1)*3-1 fem.IX(:,1)*3 fem.IX(:,2)*3-2 fem.IX(:,2)*3-1 fem.IX(:,2)*3 ...  
                fem.IX(:,3)*3-2 fem.IX(:,3)*3-1 fem.IX(:,3)*3 fem.IX(:,4)*3-2 fem.IX(:,4)*3-1 fem.IX(:,4)*3];                                    
    fem.is=reshape((kron(fem.edof,ones(1,12)))',144*fem.ne,1)';                                    
    fem.js=reshape((kron(fem.edof,ones(12,1)))',144*fem.ne,1)';    
    fem.E0 = inputs.mprop(1); fem.Emin = fem.E0*1e-9;
    fem.D = 1/((1+inputs.mprop(2))*(1-2*inputs.mprop(2))) * ...
          [1-inputs.mprop(2), inputs.mprop(2), inputs.mprop(2),  0, 0, 0;...
           inputs.mprop(2),   1-inputs.mprop(2), inputs.mprop(2),0, 0, 0;...
           inputs.mprop(2), inputs.mprop(2),  1-inputs.mprop(2), 0, 0, 0;
           0,    0,    0,    (1-2*inputs.mprop(2))/2,             0, 0;...
           0,    0,    0,    0,   (1-2*inputs.mprop(2))/2,           0;...
           0,    0,    0,    0,   0,          (1-2*inputs.mprop(2))/2];       
    nx=[fem.X(fem.IX(:,1),1) fem.X(fem.IX(:,2),1) ...
            fem.X(fem.IX(:,3),1) fem.X(fem.IX(:,4),1)]; % element node x location
    ny=[fem.X(fem.IX(:,1),2) fem.X(fem.IX(:,2),2) ...
            fem.X(fem.IX(:,3),2) fem.X(fem.IX(:,4),2)]; % element node y location   
    nz=[fem.X(fem.IX(:,1),3) fem.X(fem.IX(:,2),3) ...
            fem.X(fem.IX(:,3),3) fem.X(fem.IX(:,4),3)]; % element node z location   
    
% Build solid stiffnes matrix (K_S) and volume vector ()
    for e=1:fem.ne
        px=[nx(e,1) nx(e,2) nx(e,3) nx(e,4)]'; 
        py=[ny(e,1) ny(e,2) ny(e,3) ny(e,4)]';
        pz=[nz(e,1) nz(e,2) nz(e,3) nz(e,4)]';        
        Ve = det([ones(4,1), px, py, pz])/6;
        B = zeros(6,12);
        for i=1:4
            ind = [1:4]; ind(i) = []; % Remove i-th index
            pm =[1 -1 1 -1];
            beta  = -pm(i)*det([ones(3,1), py(ind), pz(ind)]);
            gamma = pm(i)*det([ones(3,1), px(ind), pz(ind)]);
            delta = -pm(i)*det([ones(3,1), px(ind), py(ind)]);
            B(:, (i-1)*3+1:(i-1)*3+3) = 1/(6*Ve) * [beta,  0,     0;...
                                                  0,     gamma, 0;... 
                                                  0,     0,     delta;...
                                                  gamma, beta,  0;...
                                                  0,     delta, gamma;...
                                                  delta, 0,     beta];
        end
        fem.K_S((e-1)*144+1:(e)*144) = reshape(B'*fem.D*B*Ve,144,1);
        fem.Ve(e)=Ve;
    end       
       
% Setting of FEM boundary condition variables
    for i=1:length(msh.bc.r.u); msh.bc.r.u{i}=node_ind_inv(msh.bc.r.u{i}); end
    for i=1:length(msh.bc.r.v); msh.bc.r.v{i}=node_ind_inv(msh.bc.r.v{i}); end
    for i=1:length(msh.bc.r.w); msh.bc.r.w{i}=node_ind_inv(msh.bc.r.w{i}); end
    for i=1:length(msh.bc.f); msh.bc.f{i}=node_ind_inv(msh.bc.f{i}); end
    for i=1:length(msh.lc.t); msh.lc.t{i}=node_ind_inv(msh.lc.t{i}); end
    for i=1:length(msh.lc.b); msh.lc.b{i}=node_ind_inv(msh.lc.b{i}); end       
       
    fem.bcdof = [];     fem.bcval = [];
    for i=1:length(msh.bc.r.u)  % roller boundary condition (u = 0)
        fem.bcdof = [fem.bcdof; reshape(msh.bc.r.u{i},[],1)*3-2];
        fem.bcval = [fem.bcval; zeros(numel(msh.bc.r.u{i}),1)];
    end
    for i=1:length(msh.bc.r.v)  % roller boundary condition (v = 0)
        fem.bcdof = [fem.bcdof; reshape(msh.bc.r.v{i},[],1)*3-1];
        fem.bcval = [fem.bcval; zeros(numel(msh.bc.r.v{i}),1)];
    end        
    for i=1:length(msh.bc.r.w)  % roller boundary condition (w = 0)
        fem.bcdof = [fem.bcdof; reshape(msh.bc.r.w{i},[],1)*3];
        fem.bcval = [fem.bcval; zeros(numel(msh.bc.r.w{i}),1)];
    end
    for i=1:length(msh.bc.f)   % fixed boundary condition (u=v=w = 0)
        fem.bcdof = [fem.bcdof; reshape(msh.bc.f{i},[],1)*3-2; ...
          reshape(msh.bc.f{i},[],1)*3-1; reshape(msh.bc.f{i},[],1)*3];
        fem.bcval = [fem.bcval; zeros(3*numel(msh.bc.f{i}),1)];
    end    
    
% Setting of FEM traction force variables
    fem.fdof = [];      fem.fval = [];
    for i=1:length(msh.lc.t)
        for j=1:size(msh.lc.t{i},1)
            for k=1:3
                vec_a(k)=fem.X(msh.lc.t{i}(j,3),k)-fem.X(msh.lc.t{i}(j,1),k);
                vec_b(k)=fem.X(msh.lc.t{i}(j,2),k)-fem.X(msh.lc.t{i}(j,1),k);                
            end
            S=norm(cross(vec_a,vec_b))/2;
            for k=1:3
                fem.fdof = [fem.fdof; msh.lc.t{i}(j,k)*3-2; msh.lc.t{i}(j,k)*3-1; msh.lc.t{i}(j,k)*3];
                fem.fval = [fem.fval; S/3*inputs.force(1); S/3*inputs.force(2); S/3*inputs.force(3)];
            end
        end
    end
    
end