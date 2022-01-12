function [] = fun_post(fem,opt,inputs)

    save(strcat('Result_',inputs.modelname,'.mat'));

    % Plot
    figure(4); clf(4); h=figure(4);
    set(h, 'Position', [5, 485, 1500, 800]); hold on;
    
    % Calculate fdv values at grid includinge external domain
    box=[min(fem.extX);max(fem.extX)];
    [Xgd,Ygd,Zgd]=ndgrid([box(1,1):inputs.resol:box(2,1)],...
                         [box(1,2):inputs.resol:box(2,2)],...
                         [box(1,3):inputs.resol:box(2,3)]);    
    Fi=scatteredInterpolant([fem.X;fem.extX],[opt.fdv;-1*ones(size(fem.extX,1),1)],'natural');
    Rlfdy=Fi(Xgd,Ygd,Zgd);
    
    if inputs.xmir==1
        Xgd=cat(1,flip(-Xgd,1),Xgd);
        Ygd=cat(1,Ygd,Ygd);
        Zgd=cat(1,Zgd,Zgd);
        Rlfdy=cat(1,flip(Rlfdy,1),Rlfdy);
    end
    if inputs.ymir==1    
        Xgd=cat(2,Xgd,Xgd);
        Ygd=cat(2,flip(-Ygd,2),Ygd);
        Zgd=cat(2,Zgd,Zgd);    
        Rlfdy=cat(2,flip(Rlfdy,2),Rlfdy);
    end
    if inputs.zmir==1    
        Xgd=cat(3,Xgd,Xgd);
        Ygd=cat(3,Ygd,Ygd);            
        Zgd=cat(3,flip(-Zgd,3),Zgd);
        Rlfdy=cat(3,flip(Rlfdy,3),Rlfdy);
    end
    
    stl_fv = isosurface(Xgd,Ygd,Zgd,Rlfdy,0);    
    patch(stl_fv,'FaceColor',[0.5 0.5 0.5],'Edgecolor','none');axis equal; 
    xlim([(-(inputs.xmir==1)*box(2,1)+(inputs.xmir==0)*box(1,1))*1.1 box(2,1)*1.1]);
    ylim([(-(inputs.ymir==1)*box(2,2)+(inputs.ymir==0)*box(1,2))*1.1 box(2,2)*1.1]);
    zlim([(-(inputs.zmir==1)*box(2,3)+(inputs.zmir==0)*box(1,3))*1.1 box(2,3)*1.1]);
    grid on;view([-0.5,-3.5,3]);camlight('headlight');camlight('left')
    
    stlwrite(strcat(inputs.modelname,'.stl'),stl_fv);
    
end
