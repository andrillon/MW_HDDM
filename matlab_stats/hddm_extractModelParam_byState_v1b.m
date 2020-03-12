%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))

%% Initialize variables.
for nTask=1:2
    if nTask==1
        filename = '/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Digits/Model_2/model_digit_2_stats.csv';
    else
        filename = '/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Faces/Model_2/model_face_2_stats.csv';
    end
    delimiter = ',';
    startRow = 2;
    formatSpec = '%s%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    
    model_D = table;
    model_D.VarName1 = cellstr(dataArray{:, 1});
    model_D.mean = dataArray{:, 2};
    model_D.std = dataArray{:, 3};
    model_D.q = dataArray{:, 4};
    model_D.q1 = dataArray{:, 5};
    model_D.q2 = dataArray{:, 6};
    model_D.q3 = dataArray{:, 7};
    model_D.q4 = dataArray{:, 8};
    model_D.mcerr = dataArray{:, 9};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    
     if nTask==1
     filename = '/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Digits/Model_2/model_digit_2_traces.csv';
    else
     filename = '/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Faces/Model_2/model_faces_2_traces.csv';
     end
    delimiter = ',';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    
    traces_D = table;
    traces_D.VarName1 = dataArray{:, 1};
    traces_D.a1 = dataArray{:, 2};
    traces_D.a2 = dataArray{:, 3};
    traces_D.a3 = dataArray{:, 4};
    traces_D.a_std = dataArray{:, 5};
    traces_D.t1 = dataArray{:, 6};
    traces_D.t2 = dataArray{:, 7};
    traces_D.t3 = dataArray{:, 8};
    traces_D.t_std = dataArray{:, 9};
    traces_D.v10 = dataArray{:, 10};
    traces_D.v11 = dataArray{:, 11};
    traces_D.v20 = dataArray{:, 12};
    traces_D.v21 = dataArray{:, 13};
    traces_D.v30 = dataArray{:, 14};
    traces_D.v31 = dataArray{:, 15};
    traces_D.v_std = dataArray{:, 16};
    traces_D.z_std = dataArray{:, 17};
    traces_D.z_trans1 = dataArray{:, 18};
    traces_D.z_trans2 = dataArray{:, 19};
    traces_D.z_trans3 = dataArray{:, 20};
    traces_D.z1 = dataArray{:, 21};
    traces_D.z2 = dataArray{:, 22};
    traces_D.z3 = dataArray{:, 23};
    traces_D.v_bias1 = dataArray{:, 24};
    traces_D.v_bias2 = dataArray{:, 25};
    traces_D.v_bias3 = dataArray{:, 26};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;

    %%
    var_2factor={'v'};
    var_1factor={'a','t','z'};
    
    mat_LME=[];
    for n=1:length(var_2factor)
        rows=find(~cellfun(@isempty,regexp(model_D.VarName1,sprintf('^%s_subj*',var_2factor{n}))));
        rownames=strvcat(model_D.VarName1);
        mat_LME{n}(:,1)=model_D.mean(rows);
        mat_LME{n}(:,2)=str2num(rownames(rows,8));
        mat_LME{n}(:,3)=str2num(rownames(rows,10));
        mat_LME{n}(:,4)=str2num(rownames(rows,13:15));
    end
    tbl_v=array2table(mat_LME{1},'VariableNames',{'v','state','cond','subj'});
    tbl_v.cond=categorical(tbl_v.cond);
    tbl_v.state=categorical(tbl_v.state);
    tbl_v.subj=categorical(tbl_v.subj);
    % tbl_v.state(tbl_v.state=='1')='ON';
    % tbl_v.state(tbl_v.state=='2')='MW';
    % tbl_v.state(tbl_v.state=='3')='MB';
    % tbl_v.cond(tbl_v.cond=='1')='GO';
    % tbl_v.cond(tbl_v.cond=='0')='NG';
    
    mdl.D.v.M0= fitlme(tbl_v,'v~1+cond+(1|subj)');
    mdl.D.v.M1= fitlme(tbl_v,sprintf('v~1+cond+state+(1|subj)'));
    mdl.D.v.M2= fitlme(tbl_v,sprintf('v~1+cond*state+(1|subj)'));
    
    %% rain cloud plot
    state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
    cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];
    
    rcp_data=[];
    rcp_traces=[];
    for i = 1:3
        for j = 1:2
            rcp_data{i, j} = mat_LME{1}(mat_LME{1}(:,2) == i & mat_LME{1}(:, 3) ==j-1,1);
            eval(sprintf('rcp_traces{i, j}=traces_D.v%g%g;',i,j-1));
        end
    end
    figure; set(gcf,'position',[93   310   907   488])
    for ncond=1:2
        subplot(1, 2, ncond); format_fig;
        for nstate=1:3
            h = raincloud_plot_and_posterior({rcp_data{nstate,ncond},rcp_traces{nstate,ncond}}, 'band_width',0.25, 'box_on', 1, 'color', state_colours(nstate,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', 0.4*(nstate-1)+0.1, 'dot_dodge_amount',  0.4*(nstate-1)+0.3,...
                'box_col_match', 1);
            set(h{2},'SizeData',144)
            for nstate2=nstate:3
                [pV_2factor{ncond}(nstate,nstate2)]=ranksum(rcp_data{nstate,ncond},rcp_data{nstate2,ncond});
                [pVt_2factor{ncond}(nstate,nstate2)]=ranksum(rcp_traces{nstate,ncond},rcp_traces{nstate2,ncond});
            end
        end
        % legend([h1{1} h2{1} h3{1}], {'ON', 'MW', 'MB'});
        % h   = rm_raincloud(rcp_data, state_colours);
        if ncond==1
            if nTask==1
            set(gca, 'Xlim',[-4 0.5],'YLim', [-1.1 1.1]);
            else
            set(gca, 'Xlim',[-3.5 0.5],'YLim', [-1 1]*1.2);
            end
            title(['HDDM - v - NOGO']);
        else
            if nTask==1
                set(gca, 'Xlim',[0 9],'YLim', [-1 1]*0.6);
            else
                set(gca, 'Xlim',[1 8],'YLim', [-1 1]*0.8);
            end
            title(['HDDM - v - GO']);
        end
        format_fig;
    end
    
    % u_subj=unique(tbl_v.subj);
    % mat_v_raov=nan(length(u_subj),6);
    % for ns=1:length(u_subj)
    %     if sum(tbl_v.cond=='0' & tbl_v.state=='1' & tbl_v.subj==u_subj(ns))~=0
    %         mat_v_raov(ns,1)=tbl_v.v(tbl_v.cond=='0' & tbl_v.state=='1' & tbl_v.subj==u_subj(ns));
    %     end
    %     if sum(tbl_v.cond=='0' & tbl_v.state=='2' & tbl_v.subj==u_subj(ns))~=0
    %         mat_v_raov(ns,2)=tbl_v.v(tbl_v.cond=='0' & tbl_v.state=='2' & tbl_v.subj==u_subj(ns));
    %     end
    %     if sum(tbl_v.cond=='0' & tbl_v.state=='3' & tbl_v.subj==u_subj(ns))~=0
    %         mat_v_raov(ns,3)=tbl_v.v(tbl_v.cond=='0' & tbl_v.state=='3' & tbl_v.subj==u_subj(ns));
    %     end
    %
    %     if sum(tbl_v.cond=='1' & tbl_v.state=='1' & tbl_v.subj==u_subj(ns))~=0
    %         mat_v_raov(ns,4)=tbl_v.v(tbl_v.cond=='1' & tbl_v.state=='1' & tbl_v.subj==u_subj(ns));
    %     end
    %     if sum(tbl_v.cond=='1' & tbl_v.state=='2' & tbl_v.subj==u_subj(ns))~=0
    %         mat_v_raov(ns,5)=tbl_v.v(tbl_v.cond=='1' & tbl_v.state=='2' & tbl_v.subj==u_subj(ns));
    %     end
    %     if sum(tbl_v.cond=='1' & tbl_v.state=='3' & tbl_v.subj==u_subj(ns))~=0
    %         mat_v_raov(ns,6)=tbl_v.v(tbl_v.cond=='1' & tbl_v.state=='3' & tbl_v.subj==u_subj(ns));
    %     end
    % end
    % writetable(tbl_v,['/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Digits/Model_2/model_digit_2_v.csv']);
    %%
    mat_LME2=[];
    for n=1:length(var_1factor)
        rows=find(~cellfun(@isempty,regexp(model_D.VarName1,sprintf('^%s_subj*',var_1factor{n}))));
        rownames=strvcat(model_D.VarName1);
        mat_LME2{n}(:,1)=model_D.mean(rows);
        mat_LME2{n}(:,2)=str2num(rownames(rows,8));
        mat_LME2{n}(:,3)=str2num(rownames(rows,11:13));
    end
    tbl_a=array2table(mat_LME2{1},'VariableNames',{'a','state','subj'});
    tbl_a.state=categorical(tbl_a.state);
    tbl_a.subj=categorical(tbl_a.subj);
    
    tbl_t=array2table(mat_LME2{2},'VariableNames',{'t','state','subj'});
    tbl_t.state=categorical(tbl_t.state);
    tbl_t.subj=categorical(tbl_t.subj);
    
    tbl_z=array2table(mat_LME2{3},'VariableNames',{'z','state','subj'});
    tbl_z.state=categorical(tbl_z.state);
    tbl_z.subj=categorical(tbl_z.subj);
    
    mdl.D.a.M0= fitlme(tbl_a,'a~1+(1|subj)');
    mdl.D.a.M1= fitlme(tbl_a,sprintf('a~1+state+(1|subj)'));
    
    mdl.D.t.M0= fitlme(tbl_t,'t~1+(1|subj)');
    mdl.D.t.M1= fitlme(tbl_t,sprintf('t~1+state+(1|subj)'));
    
    mdl.D.z.M0= fitlme(tbl_z,'z~1+(1|subj)');
    mdl.D.z.M1= fitlme(tbl_z,sprintf('z~1+state+(1|subj)'));
    
    %%
    rcp_data=[];
    rcp_traces=[];
    for i = 1:3
        for j=1:3
            rcp_data{i, j} = mat_LME2{j}(mat_LME2{j}(:,2) == i,1);
            eval(sprintf('rcp_traces{i, j}=traces_D.%s%g;',var_1factor{j},i));
                        
        end
    end
    figure; set(gcf,'position',[93   310   1360   488])
    for ncond=1:3
        subplot(1, 3, ncond); format_fig;
        for nstate=1:3
            h = raincloud_plot_and_posterior({rcp_data{nstate,ncond},rcp_traces{nstate,ncond}}, 'box_on', 1, 'color', state_colours(nstate,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', 0.4*(nstate-1)+0.1, 'dot_dodge_amount',  0.4*(nstate-1)+0.3,...
                'box_col_match', 1);
            set(h{2},'SizeData',144)
            for nstate2=nstate:3
                [pV_1factor{ncond}(nstate,nstate2)]=ranksum(rcp_data{nstate,ncond},rcp_data{nstate2,ncond});
                [pVt_1factor{ncond}(nstate,nstate2)]=ranksum(rcp_traces{nstate,ncond},rcp_traces{nstate2,ncond});
            end
        end
        % legend([h1{1} h2{1} h3{1}], {'ON', 'MW', 'MB'});
        % h   = rm_raincloud(rcp_data, state_colours);
        if ncond==1
            if nTask==1
                set(gca, 'XLim', [0.2 3], 'YLim', [-1 1]*2);
            else
                set(gca, 'XLim', [1.2 3.2], 'YLim', [-1 1]*8);
            end
            title(['HDDM - a']);
        elseif ncond==2
            if nTask==1
                set(gca, 'XLim', [0.15 .4], 'YLim', [-1 1]*20);
            else
                set(gca, 'XLim', [0.1 0.35], 'YLim', [-1 1]*50);
            end
            title(['HDDM - t']);
        elseif ncond==3
            if nTask==1
                set(gca, 'XLim', [0.15 0.3], 'YLim', [-1 1]*10);
            else
                set(gca, 'XLim', [0.1 0.35], 'YLim', [-1 1]*50);
            end
            title(['HDDM - z']);
        end
        format_fig;
    end
end