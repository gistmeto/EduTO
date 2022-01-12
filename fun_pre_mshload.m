function [msh] = fun_pre_mshload(modelname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective
%   Convert "filename.msh(v4.1)" to "Inner FEM data(IX,X,bc,lc)"
%
% Requirement
%   (GMSH) Physical tag is required to identify domain(Design & Non-design)
%           and condition(boundary & loading)
%
%  ***Domain identification***
%
%  If there are no non-design domain, any physical tag is not required to 
%  identify design domain. Because default domain setting is design domain.
%  However if there are a non-design domian, physical tag is required. A
%  physical tag must contain specific keyword for target domain.
%
%  'Design', 'NonDesign', and 'External' are keywords for check non-design 
%  domain.                               
%  Ex) 'Design2', 'NonDesign1', 'External' ...
%  ***************************                                             
%
%
%  ***Condition identification***
%
%  Physical tag are required to identify conditions. A physical tag must
%  contain specific keyword for target condition.
%
%  'RollerU','RollerV','RollerW', 'Fixed' are keywords for check boundary 
%      condition
%  Ex) RollerU1, RollerW2, ...
%      Fixed1, Fixed2, ...
%
%  'Traction', 'Body' are keywords for check loading condition             
%  Ex) 'Traction1', 'Body3', 'Body3' ... 
%  ******************************
%  
% IX : element conectivity + Domain index information (tetrahedron only)
%  Type : index matrix
%    IX(:,1) : 1st node index
%    IX(:,2) : 2nd node index
%    IX(:,3) : 3rd node index
%    IX(:,4) : 4th node index
%    IX(:,5) : Domain index      % 1 : Design // 2 : Nondesign // 3 : External
%  
% X : Node coordinates information
%  Type : coordinates matrix
%    X(:,1) : X-axis coordinated
%    X(:,2) : Y-axis coordinated
%    X(:,3) : Z-axis coordinated
%  
% bc : Bondary condition
%  Type : structure
%    bc.r : Roller condition                                               
%      Type : {Cell}(elemenet conectivity matrix)
%        bc.r.u{ind} : if roller condition is u=0
%        bc.r.v{ind} : if roller condition is v=0
%        bc.r.w{ind} : if roller condition is w=0
%    bc.f : Fixed condition
%      Type : index vector
%   
% lc : Loading condition {Cell}(elemenet conectivity matrix)
%   lc.t{ind} : Traction loading
%   lc.b{ind} : Bady loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fid = fopen([modelname,'.msh']);
flag_category = 0; % Flag for ctaegroy
flag_state = 0;    % Flag for state

while 1
    % File reading part
    tline = fgetl(fid);
    if ~ischar(tline), break, end % if end the file, break the while loop
    
    % line 
    [data,~,errmsg,nextindex] = sscanf(tline,'%f');
%     disp(tline)
%     disp(data)
    
    if (~isempty(errmsg)) && (nextindex == 1)
        if strcmp(tline,'$PhysicalNames')
          flag_category = 1;
        elseif strcmp(tline,'$Entities')
          flag_category = 2;
        elseif strcmp(tline,'$Nodes')
          flag_category = 4;
        elseif strcmp(tline,'$Elements')
          flag_category = 5;
        end
        flag_state = 0;
        continue
    end
    
    if flag_category == 1 % $PhysicalNames
        if flag_state == 0
            if length(data) ~= 1
                error
            end
            num_phys_entities = data;
            flag_state = 1;
            phys_tag_str.RollerU = [];                                     % TODO : Arbitrary roller surface implementation
            phys_tag_str.RollerV = [];
            phys_tag_str.RollerW = [];
            phys_tag_str.Fixed = [];
            phys_tag_str.Traction = [];
            phys_tag_str.Body = [];
            phys_tag_str.Design = [];
            phys_tag_str.NonDesign = [];
            phys_tag_str.External = [];
            ind_phys_entities = 0;
        
                                                                           % TODO : check Matlab version dependency
        elseif flag_state == 1  % TODO : Check error
            if (strfind(tline(nextindex:end),'"RollerU"') == 1) 
                ind_phys_entities = ind_phys_entities + 1;
                Boundary_roller_u_list = [];
                phys_tag_str.RollerU = [phys_tag_str.RollerU, data(2)];
            elseif (strfind(tline(nextindex:end),'"RollerV"') == 1) 
                ind_phys_entities = ind_phys_entities + 1;
                Boundary_roller_v_list = [];
                phys_tag_str.RollerV = [phys_tag_str.RollerV, data(2)];
            elseif (strfind(tline(nextindex:end),'"RollerW"') == 1) 
                ind_phys_entities = ind_phys_entities + 1;
                Boundary_roller_w_list = [];
                phys_tag_str.RollerW = [phys_tag_str.RollerW, data(2)];
            elseif (strfind(tline(nextindex:end),'"Fixed"') == 1)
                ind_phys_entities = ind_phys_entities + 1;
                Boundary_fixed_list = [];
                phys_tag_str.Fixed = [phys_tag_str.Fixed, data(2)];
            elseif (strfind(tline(nextindex:end),'"Traction"') == 1)
                Traction_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Traction = [phys_tag_str.Traction, data(2)];
            elseif (strfind(tline(nextindex:end),'"Body"') == 1)
                Body_list = [];                
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Body = [phys_tag_str.Body, data(2)];
            elseif (strfind(tline(nextindex:end),'"Design"') == 1)
                Design_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Design = [phys_tag_str.Design, data(2)];
            elseif ((strfind(tline(nextindex:end),'"NonDesign"')) == 1)
                Nondesign_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.NonDesign = [phys_tag_str.NonDesign, data(2)];
            elseif (strfind(tline(nextindex:end),'"External"') == 1)
                External_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.External = [phys_tag_str.External, data(2)];
            end
        end
        
        if ind_phys_entities >= num_phys_entities
            flag_state = 0;
            ind_phys_entities = 0;
        end
                    
    elseif flag_category == 2 % $Entities
        if flag_state == 0
            num_point = data(1);
            num_line = data(2);
            num_surf = data(3);
            num_vol = data(4);
            if exist('phys_tag_str','var')
                if ~(isempty(phys_tag_str))
                    flag_state = 1;
                    flag_entity = 0;
                    ind_entity = 0;
                end
            end
            
        elseif flag_state == 1
            if flag_entity == 0 % Point
                ind_entity = ind_entity + 1;
                if ind_entity >= num_point
                    flag_entity = 1;
                    ind_entity = 0;
                end
            elseif flag_entity == 1 % Curve
                ind_entity = ind_entity + 1;
                if ind_entity >= num_line
                    flag_entity = 2;
                    ind_entity = 0;
                end
            elseif flag_entity == 2 % Surface
                ind_entity = ind_entity + 1;
                
                if data(8) >= 1
                    for i = 1:data(8)
                        if (sum(data(8+i) == phys_tag_str.RollerU) >= 1)     % Roller U
                            Boundary_roller_u_list = [Boundary_roller_u_list; data(1)];
                        elseif (sum(data(8+i) == phys_tag_str.RollerV) >= 1)     % Roller V
                            Boundary_roller_v_list = [Boundary_roller_v_list; data(1)];
                        elseif (sum(data(8+i) == phys_tag_str.RollerW) >= 1)     % Roller W
                            Boundary_roller_w_list = [Boundary_roller_w_list; data(1)];
                        elseif (sum(data(8+i) == phys_tag_str.Fixed) >= 1)     % Fixed
                            Boundary_fixed_list = [Boundary_fixed_list; data(1)];
                        elseif (sum(data(8+i) == phys_tag_str.Traction) >= 1) % Traction
                            Traction_list = [Traction_list; data(1)];
                        end
                    end
                end
                
                if ind_entity >= num_surf
                    flag_entity = 3;
                    ind_entity = 0;
                end
            elseif flag_entity == 3 % Volume
                ind_entity = ind_entity + 1;
                
                for i = 1:data(8)
                    if (sum(data(8+i) == phys_tag_str.Body) >= 1)     % Body
                        Body_list = [Body_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.Design) >= 1)     % Design
                        Design_list = [Design_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.NonDesign) >= 1)     % Nondesign
                        Nondesign_list = [Nondesign_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.External) >= 1)     % External
                        External_list = [External_list; data(1)];
                    end
                end
                
                if ind_entity >= num_vol
                    flag_entity = 4;
                    ind_entity = 0;
                end
            end
        end
        
    elseif flag_category == 4 % $Nodes
        if flag_state == 0 % Initial state
            num_node = data(2);    % Total number of node
            p = zeros(num_node,3);  % Node coordinates data
            flag_state = 1; %
        
        elseif flag_state == 1  % Info data parsing
            num_loop = data(4);    % Number of node on this block
            nodeTag_list = zeros(num_loop,1);
            if num_loop == 0
                continue
            else
                flag_state = 2;
                ind_inner_loop = 0;
            end
        
        elseif flag_state == 2  % ind data parsing
            ind_inner_loop = ind_inner_loop + 1;
            nodeTag_list(ind_inner_loop) = data;
            if ind_inner_loop >= num_loop % If ind data parsing is ended
                flag_state = 3;
                ind_inner_loop = 0;
            end
            
        elseif flag_state == 3  % point data parsing
            ind_inner_loop = ind_inner_loop + 1;
            try
                p(nodeTag_list(ind_inner_loop),:) = data;
            catch
                disp(nodeTag_list(ind_inner_loop))
            end
            
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
            end
        elseif flag_state == 99 % error exception
            ind_inner_loop = ind_inner_loop + 1;
            
            if ind_inner_loop >= num_loop*2 % If ind data parsing is ended
                flag_state = 1;
                ind_inner_loop = 0;
            end
        end
    elseif flag_category == 5 % $Elements
        if flag_state == 0 % Initial state
            if length(data) ~= 4
                error
            end
            flag_state = 1;     %
            t = [];             % Element conectivity information
            
            bc.r.u = {};        % Boundary condition u
            bcru_ind = 0;
            bc.r.v = {};        % Boundary condition v
            bcrv_ind = 0;
            bc.r.w = {};        % Boundary condition w
            bcrw_ind = 0;
            bc.f = {};          % Boundary condition fixed
            bcf_ind = 0;
            
            lc.t = {};          % Tranction force element conectivity infromation
            lct_ind = 0;
            lc.b = {};          % Body force element conectivity information
            lcb_ind = 0;
            
        elseif flag_state == 1  % Info data parsing
            num_loop = data(4);    % Number of node on this block
            flag_state = 99;
            ind_inner_loop = 0;
            
            flag_condition = zeros(6,1); % index : 1:u,2:v,3:w,4:f,5:t,6:b
            
            if data(1) == 1      % 1D element : Point
                continue
            elseif data(1) == 2   % 2D element
                if data(3) == 2  % 1D Line
                    if exist('Boundary_roller_u_list','var')               % TODO : t_temp split
                        if (sum(Boundary_roller_u_list==data(2)) >= 1) % If data(2) is in Boundary_roller_u_list
                           flag_state = 2;
                           flag_condition(1) = 1;
                           bcru_ind = bcru_ind + 1;
                           num_elem = data(3);
                           t_temp = zeros(num_elem,3);
                        end
                    end
                    if exist('Boundary_roller_v_list','var')
                        if (sum(Boundary_roller_v_list==data(2)) >= 1)
                           flag_state = 2;
                           flag_condition(2) = 1;
                           bcrv_ind = bcrv_ind + 1;
                           num_elem = data(3);
                           t_temp = zeros(num_elem,3);
                        end
                    end
                    if exist('Boundary_roller_w_list','var')
                        if (sum(Boundary_roller_w_list==data(2)) >= 1)
                           flag_state = 2;
                           flag_condition(3) = 1;
                           bcrw_ind = bcrw_ind + 1;
                           num_elem = data(3);
                           t_temp = zeros(num_elem,3);
                        end
                    end    
                    if exist('Boundary_fixed_list','var')
                        if (sum(Boundary_fixed_list==data(2)) >= 1)
                           flag_state = 2;
                           flag_condition(4) = 1;
                           bcf_ind = bcf_ind + 1;
                           num_elem = data(3);
                           t_temp = zeros(num_elem,3);
                        end
                    end
                    if exist('Traction_list','var')
                        if (sum(Traction_list==data(2)) >= 1)
                           flag_state = 3;
                           flag_condition(5) = 1;
                           lct_ind = lct_ind+1;
                           num_elem = data(3);
                           t_temp = zeros(num_elem,3);
                        end
                    end
                    
                elseif data(3) == 3 % 2D quadrangle element
                    continue
                end
                
            elseif data(1) == 3   % 3D element
                if data(3) == 4   % 3D tetrahedron element
                    flag_state = 4;
                    num_elem = data(4);
                    t_temp = zeros(num_elem,5);
                    
                    if exist('Body_list','var')
                        if sum(Body_list==data(2))
                            flag_condition(6) = 1;
                            lcb_ind = lcb_ind+1;
                        end
                    end
                    
                    if exist('Design_list','var')
                        if sum(Design_list==data(2))
                            t_temp(:,5) = ones(num_elem,1);                % TODO : change index
                        end
                    end
                    
                    if exist('Nondesign_list','var')
                        if sum(Nondesign_list==data(2))
                            t_temp(:,5) = ones(num_elem,1)*2;
                        end
                    end
                    
                    if exist('External_list','var')
                        if sum(External_list==data(2))
                            t_temp(:,5) = ones(num_elem,1)*3;
                        end
                    end
                    
                elseif data(3) == 5 % 3D hexahedron element
                    continue
                end
            end
            
        elseif flag_state == 2  % 2D boundary element conectivity data parsing
            ind_inner_loop = ind_inner_loop + 1;
            
            t_temp(ind_inner_loop,1:3) = data(2:4);
            
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
                
                if     flag_condition(1) == 1       % u
                    bc.r.u{bcru_ind} = t_temp;
                elseif flag_condition(2) == 1       % v
                    bc.r.v{bcrv_ind} = t_temp;
                elseif flag_condition(3) == 1       % w
                    bc.r.w{bcrw_ind} = t_temp;
                elseif flag_condition(4) == 1       % f
                    bc.f{bcf_ind} = t_temp;
                end
            end
            
        elseif flag_state == 3  % 2D triangle element conectivity data parsing
            ind_inner_loop = ind_inner_loop + 1;
            t_temp(ind_inner_loop,1:3) = data(2:4);
            
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
                
                if flag_condition(5) == 1           % 
                    lc.t{lct_ind} = t_temp;
                end
            end
            
        elseif flag_state == 4  % 3D tetrahedron element conectivity data parsing
            ind_inner_loop = ind_inner_loop + 1;
            t_temp(ind_inner_loop,1:4) = data(2:5);
            
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
                t = [t; t_temp];
                
                if flag_condition(6) == 1           %
                   lc.b{lcb_ind}  = t_temp(:,1:4);
                end
            end
            
        elseif flag_state == 99 % error exception
            ind_inner_loop = ind_inner_loop + 1;
            
            if ind_inner_loop >= num_loop % If ind data parsing is ended
                flag_state = 1;
                ind_inner_loop = 0;
            end
        end
    end
end
fclose(fid);

%% Check X and IX
temp = unique(t);

if max(temp(2:end)-temp(1:end-1)) > 1
    error('X order has an error')
end

p = p(temp,:); % Remove useless node

%% Output setting
msh.IX = t;
msh.X = p;

msh.bc=bc;
msh.lc=lc;

%% Post processing
while 0
    % Plotting mesh - tetramesh
    figure(1);
    tetramesh(fem.msh_IX(:,1:4),fem.msh_X,'EdgeColor','k','FaceColor','k','FaceAlpha',0.05)
    view([-0.5,-3,3.5]); axis equal; grid on;drawnow;

    % Check design and non-deisgn domain
    figure(2); hold on;
    Design_ind = find(fem.msh_IX(:,5) == 1);
    tetramesh(fem.msh_IX(Design_ind,1:4),fem.msh_X,'EdgeColor','k','FaceColor','r','FaceAlpha',0.05)
    Nondesign_ind = find(fem.msh_IX(:,5) == 2);
    tetramesh(fem.msh_IX(Nondesign_ind,1:4),fem.msh_X,'EdgeColor','k','FaceColor','b','FaceAlpha',0.5)
    External_ind = find(fem.msh_IX(:,5) == 3);
    tetramesh(fem.msh_IX(External_ind,1:4),fem.msh_X,'EdgeColor','k','FaceColor','g','FaceAlpha',0.5)
    view([-0.5,-3,3.5]); axis equal; grid on;drawnow;    

    % Check boundary condition
    figure(3); hold on;
    tetramesh(fem.msh_IX(:,1:4),fem.msh_X,'EdgeColor','k','FaceColor','k','FaceAlpha',0.05) % Base mesh
    for i = 1:length(bc.r.u)
        trimesh(bc.r.u{i},fem.msh_X(:,1),fem.msh_X(:,2),fem.msh_X(:,3),'FaceColor','r') % Roller U
    end
    for i = 1:length(bc.r.v)
        trimesh(bc.r.v{i},fem.msh_X(:,1),fem.msh_X(:,2),fem.msh_X(:,3),'FaceColor','g') % Roller V
    end
    for i = 1:length(bc.r.w)
        trimesh(bc.r.w{i},fem.msh_X(:,1),fem.msh_X(:,2),fem.msh_X(:,3),'FaceColor','b') % Roller W
    end
    for i = 1:length(bc.f)
        trimesh(bc.f{i},fem.msh_X(:,1),fem.msh_X(:,2),fem.msh_X(:,3),'FaceColor','m') % Fixed
    end
    view([-0.5,-3,3.5]); axis equal; grid on;drawnow;    

    % Check loading condition
    figure(4); hold on;
    tetramesh(fem.msh_IX(:,1:4),fem.msh_X,'EdgeColor','k','FaceColor','k','FaceAlpha',0) % Base mesh
    for i = 1:length(lc.t)
        trimesh(lc.t{i},fem.msh_X(:,1),fem.msh_X(:,2),fem.msh_X(:,3),'FaceColor','b') % Taction loading condition
    end
    for i = 1:length(lc.b)
        tetramesh(lc.b{i},fem.msh_X,'EdgeColor','k','FaceColor','r','FaceAlpha',0.05) % Body loading condition
    end
    view([-0.5,-3,3.5]); axis equal; grid on; drawnow;
    
    break;
end
end