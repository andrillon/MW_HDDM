%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))

%% Initialize variables.
filename = '/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Omni/Model_4/model_omni_4_stats.csv';

delimiter = ',';
startRow = 2;
formatSpec = '%s%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

model_B = table;
model_B.VarName1 = cellstr(dataArray{:, 1});
model_B.mean = dataArray{:, 2};
model_B.std = dataArray{:, 3};
model_B.q = dataArray{:, 4};
model_B.q1 = dataArray{:, 5};
model_B.q2 = dataArray{:, 6};
model_B.q3 = dataArray{:, 7};
model_B.q4 = dataArray{:, 8};
model_B.mcerr = dataArray{:, 9};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;




%%
var_2factor={'v'};

mat_LME=[];
for n=1:length(var_2factor)
    rows=find(~cellfun(@isempty,regexp(model_B.VarName1,sprintf('^%s_subj*',var_2factor{n}))));
    rownames=strvcat(model_B.VarName1);
    mat_LME{n}(:,1)=model_B.mean(rows);
    mat_LME{n}(:,2)=str2num(rownames(rows,8));
    mat_LME{n}(:,3)=str2num(rownames(rows,10));
    mat_LME{n}(:,4)=str2num(rownames(rows,12));
    mat_LME{n}(:,5)=str2num(rownames(rows,15:17));
end
tbl_v=array2table(mat_LME{1},'VariableNames',{'v','state','task','cond','subj'});
tbl_v.cond=categorical(tbl_v.cond);
tbl_v.state=categorical(tbl_v.state);
tbl_v.subj=categorical(tbl_v.subj);
tbl_v.task=categorical(tbl_v.task);
% tbl_v.state(tbl_v.state=='1')='ON';
% tbl_v.state(tbl_v.state=='2')='MW';
% tbl_v.state(tbl_v.state=='3')='MB';
% tbl_v.cond(tbl_v.cond=='1')='GO';
% tbl_v.cond(tbl_v.cond=='0')='NG';

mdl.vGO.M0= fitlme(tbl_v(tbl_v.cond=='1',:),'v~1+(1|subj)');
mdl.vGO.M1= fitlme(tbl_v(tbl_v.cond=='1',:),'v~1+task+(1|subj)');
mdl.vGO.M2= fitlme(tbl_v(tbl_v.cond=='1',:),'v~1+task+state+(1|subj)'); % winning model
% mdl.vGO.M3= fitlme(tbl_v(tbl_v.cond=='1',:),'v~1+task*state+(1|subj)');

mdl.vNGO.M0= fitlme(tbl_v(tbl_v.cond=='0',:),'v~1+(1|subj)');
% mdl.vNGO.M1= fitlme(tbl_v(tbl_v.cond=='0',:),'v~1+task+(1|subj)');
mdl.vNGO.M2= fitlme(tbl_v(tbl_v.cond=='0',:),'v~1+state+(1|subj)'); % winning model
% mdl.vNGO.M3= fitlme(tbl_v(tbl_v.cond=='0',:),'v~1+task+state+(1|subj)');


% mdl.vbias.M0= fitlme(tbl_vbias,'v~1+(1|subj)');
% mdl.vbias.M1= fitlme(tbl_vbias,'v~1+task+(1|subj)');
% mdl.vbias.M2= fitlme(tbl_vbias,'v~1+task+state+(1|subj)');
% mdl.vbias.M3= fitlme(tbl_vbias,'v~1+task*state+(1|subj)');

%% rain cloud plot
state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

for ntask=1:2
    rcp_data=[];
    for i = 1:3
        for j = 1:2
            rcp_data{i, j} = mat_LME{1}(mat_LME{1}(:,2) == i & mat_LME{1}(:, 4) ==j-1 & mat_LME{1}(:, 3) ==ntask,1);
        end
    end
    figure; set(gcf,'position',[93   310   880   880])
    for ncond=1:2
        subplot(1, 2, ncond); format_fig;
        for nstate=1:3
            h = raincloud_plot(rcp_data{nstate,ncond}, 'band_width',0.25, 'box_on', 1, 'color', state_colours(nstate,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', 0.4*(nstate-1)+0.1, 'dot_dodge_amount',  0.4*(nstate-1)+0.3,...
                'box_col_match', 1);
            set(h{2},'SizeData',144)
            for nstate2=nstate:3
                [pV_2factor{ncond}(nstate,nstate2)]=ranksum(rcp_data{nstate,ncond},rcp_data{nstate2,ncond});
            end
        end
        % legend([h1{1} h2{1} h3{1}], {'ON', 'MW', 'MB'});
        % h   = rm_raincloud(rcp_data, state_colours);
        if ncond==1
            if ntask==1
                set(gca, 'Xlim',[-4 0.5],'YLim', [-1 1]*1.15);
                title(['D - v - NOGO']);
            else
                set(gca, 'Xlim',[-4 0.5],'YLim', [-1 1]*1.2);
                title(['F - v - NOGO']);
            end
        else
            if ntask==1
                set(gca, 'Xlim',[1 8],'YLim', [-1 1]*0.9);
                title(['D - v - GO']);
            else
                set(gca, 'Xlim',[1 8],'YLim', [-1 1]*0.8);
                title(['F - v - GO']);
            end
        end
        format_fig;
    end
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
var_1factor={'a','t','z'};
mat_LME2=[];
for n=1:length(var_1factor)
    rows=find(~cellfun(@isempty,regexp(model_B.VarName1,sprintf('^%s_subj*',var_1factor{n}))));
    rownames=strvcat(model_B.VarName1);
    mat_LME2{n}(:,1)=model_B.mean(rows);
    mat_LME2{n}(:,2)=str2num(rownames(rows,8));
    mat_LME2{n}(:,3)=str2num(rownames(rows,10));
    mat_LME2{n}(:,4)=str2num(rownames(rows,13:15));
end
tbl_a=array2table(mat_LME2{1},'VariableNames',{'a','state','task','subj'});
tbl_a.state=categorical(tbl_a.state);
tbl_a.subj=categorical(tbl_a.subj);

tbl_t=array2table(mat_LME2{2},'VariableNames',{'t','state','task','subj'});
tbl_t.state=categorical(tbl_t.state);
tbl_t.subj=categorical(tbl_t.subj);

tbl_z=array2table(mat_LME2{3},'VariableNames',{'z','state','task','subj'});
tbl_z.state=categorical(tbl_z.state);
tbl_z.subj=categorical(tbl_z.subj);

mdl.a.M0= fitlme(tbl_a,'a~1+(1|subj)');
mdl.a.M1= fitlme(tbl_a,sprintf('a~1+task+(1|subj)'));
mdl.a.M2= fitlme(tbl_a,sprintf('a~1+task+state+(1|subj)')); % winning model
% mdl.a.M3= fitlme(tbl_a,sprintf('a~1+task*state+(1|subj)'));

mdl.t.M0= fitlme(tbl_t,'t~1+(1|subj)');
mdl.t.M1= fitlme(tbl_t,sprintf('t~1+state+(1|subj)'));
mdl.t.M2= fitlme(tbl_t,sprintf('t~1+task+state+(1|subj)'));
mdl.t.M3= fitlme(tbl_t,sprintf('t~1+task*state+(1|subj)')); % winning model


mdl.z.M0= fitlme(tbl_z,'z~1+(1|subj)');
mdl.z.M1= fitlme(tbl_z,sprintf('z~1+state+(1|subj)'));
mdl.z.M2= fitlme(tbl_z,sprintf('z~1+task+state+(1|subj)')); % winning model
% mdl.z.M3= fitlme(tbl_z,sprintf('z~1+task*state+(1|subj)'));


%%
for ntask=1:2
    rcp_data=[];
    for i = 1:3
        for j=1:3
            rcp_data{i, j} = mat_LME2{j}(mat_LME2{j}(:,2) == i & mat_LME2{j}(:,3) == ntask,1);
        end
    end
    %     figure; set(gcf,'position',[93   310   1360   488])
        figure; set(gcf,'position',[93   310   880   880])
for ncond=1:3
        subplot(1, 3, ncond); format_fig;
        for nstate=1:3
            h = raincloud_plot(rcp_data{nstate,ncond}, 'band_width',[], 'box_on', 1, 'color', state_colours(nstate,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', 0.4*(nstate-1)+0.1, 'dot_dodge_amount',  0.4*(nstate-1)+0.3,...
                'box_col_match', 1);
            set(h{2},'SizeData',144)
            for nstate2=nstate:3
                [pV_1factor{ncond}(nstate,nstate2)]=ranksum(rcp_data{nstate,ncond},rcp_data{nstate2,ncond});
            end
        end
        % legend([h1{1} h2{1} h3{1}], {'ON', 'MW', 'MB'});
        % h   = rm_raincloud(rcp_data, state_colours);
        if ncond==1
            %             if nTask==1
            %                 set(gca, 'XLim', [0.2 3], 'YLim', [-1 1]*2);
            %             else
            %                 set(gca, 'XLim', [1 3.5], 'YLim', [-1 1]*1.5);
            %             end
            %             title(['HDDM - a']);
            %         elseif ncond==2
            %         if ntask==1
            %             set(gca, 'XLim', [0.15 .4], 'YLim', [-1 1]*20);
            %         else
            %             set(gca, 'XLim', [0.1 0.35], 'YLim', [-1 1]*16);
            %         end
            title(['HDDM - t']);
            %         elseif ncond==3
            %             if nTask==1
            %                 set(gca, 'XLim', [0.1 0.55], 'YLim', [-1 1]*10);
            %                 title(['D - z']);
            %             else
            %                 set(gca, 'XLim', [0.1 0.55], 'YLim', [-1 1]*12.5);
            %                 title(['F - z']);
            %             end
        end
        format_fig;
    end
    % if nTask==1
    %     export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Digits_perState_noZ.fig'])
    %     export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Digits_perState_noZ.eps'],'-r 300')
    % else
    %     export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Faces_perState_noZ.fig'])
    %     export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Faces_perState_noZ.eps'],'-r 300')
    %
    % end
end