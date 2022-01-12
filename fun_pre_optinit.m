function [opt] = fun_pre_optinit(inputs,fem)

% Separation of design and non-design domains
    nd_ele = find(fem.IX(:,5)==2);           
    opt.dof_nd = unique([fem.IX(nd_ele,1);fem.IX(nd_ele,2);fem.IX(nd_ele,3);fem.IX(nd_ele,4)]);
    opt.dof_dd = setdiff([1:fem.nn]',opt.dof_nd);   
    
% Setting problem parameters
    opt.VT = inputs.VT;    opt.volfrac = inputs.volfrac;
    opt.penal = inputs.penal;
    
% Setting of optimization variables    
    opt.dv = ones(length(opt.dof_dd),1) * inputs.initdv;        
    opt.nv(opt.dof_nd,1) = 1;       % Varabile including non-deisgn domain
    opt.nv(opt.dof_dd,1) = opt.dv;
    
    opt.dvold=opt.dv;           opt.dvolder=opt.dv;            
    opt.dvmin=opt.dv*0-1;       opt.dvmax=opt.dv*0+1;

    opt.iter=1;                 opt.deltaf = 1.0;
    
    opt.bt = inputs.bt;

% Setting of MMA variables
    MMA.a0 = 1;                 MMA.a = [0;];   
    MMA.c = [inputs.MMA_c;];     MMA.d = [1;];
    MMA.low = opt.dvmin;        MMA.upp = opt.dvmax;
    
    opt.MMA = MMA;

% Build filter element stiffness (Kft) and Transformation matrix (Tft)
    nx=[fem.X(fem.IX(:,1),1) fem.X(fem.IX(:,2),1) ...
            fem.X(fem.IX(:,3),1) fem.X(fem.IX(:,4),1)]; % element node x location
    ny=[fem.X(fem.IX(:,1),2) fem.X(fem.IX(:,2),2) ...
            fem.X(fem.IX(:,3),2) fem.X(fem.IX(:,4),2)]; % element node y location   
    nz=[fem.X(fem.IX(:,1),3) fem.X(fem.IX(:,2),3) ...
            fem.X(fem.IX(:,3),3) fem.X(fem.IX(:,4),3)]; % element node z location   
    Kd = (inputs.rmin/2/sqrt(3))^2*[1 0 0; 0 1 0; 0 0 1];
    NN = [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2]/12;
    isf=reshape((kron(fem.IX(:,1:4),ones(1,4)))',16*fem.ne,1)';                                    
    jsf=reshape((kron(fem.IX(:,1:4),ones(4,1)))',16*fem.ne,1)';
    for e=1:fem.ne
        px=[nx(e,1) nx(e,2) nx(e,3) nx(e,4)]'; 
        py=[ny(e,1) ny(e,2) ny(e,3) ny(e,4)]';
        pz=[nz(e,1) nz(e,2) nz(e,3) nz(e,4)]';       
        Ve = fem.Ve(e);
        B = zeros(3,4);
        for i=1:4
            ind = [1:4]; ind(i) = []; % Remove i-th index
            pm =[1 -1 1 -1];
            beta  = -pm(i)*det([ones(3,1), py(ind), pz(ind)]);
            gamma = pm(i)*det([ones(3,1), px(ind), pz(ind)]);
            delta = -pm(i)*det([ones(3,1), px(ind), py(ind)]);
            B(:,i) = 1/(6*Ve) * [beta; gamma; delta];
        end
        Kft((e-1)*16+1:(e)*16)=reshape((B'*Kd*B + NN)*Ve,16,1);
        Tft((e-1)*16+1:(e)*16)=reshape(NN*Ve,16,1);
    end
    opt.chol_Kft = chol(sparse(isf,jsf,Kft),'lower');
    opt.Tft=sparse(isf,jsf,Tft);

% Build matrix for transformation from nodal to element density (Ten)
    opt.Ten=sparse(fem.ne,fem.nn);
    for e=1:fem.ne
       opt.Ten(e,fem.IX(e,1:4))=1/4;
    end
end