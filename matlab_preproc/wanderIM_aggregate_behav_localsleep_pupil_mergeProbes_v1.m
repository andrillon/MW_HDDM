%% Init
clear all
close all

% local file, set up paths
run ../localdef_wanderIM
addpath(genpath(lscpTools_path))

data_path=[root_path filesep 'behav/'];
eeg_path2=[root_path filesep 'preproc_ica'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);
load([pwd filesep '..' filesep '..' filesep 'WanderIM' filesep 'paper' filesep 'paper_SubID'])

%%
table_res=[];
for n=1:length(files)
    % load behav
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    if ~ismember(SubID,GoodSudID) || exist([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v4.mat'])==0
        continue;
    end
   
    % load matrices
    load([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v4'])
   
    % compute prctile pupil
    Pup_F=hddm_subj(hddm_subj(:,6)==1,end);
    boundaries=prctile(Pup_F,0:20:100);
    prPup_F=nan(size(Pup_F,1),1);
    for nPer=1:5
        prPup_F(Pup_F>=boundaries(nPer) & Pup_F<=boundaries(nPer+1))=nPer;
    end
    
    Pup_D=hddm_subj(hddm_subj(:,6)==2,end);
    boundaries=prctile(Pup_D,0:20:100);
    prPup_D=nan(size(Pup_D,1),1);
    for nPer=1:5
        prPup_D(Pup_D>=boundaries(nPer) & Pup_D<=boundaries(nPer+1))=nPer;
    end
    
    new_Pup=hddm_subj(:,end);
    new_Pup(hddm_subj(:,6)==1)=prPup_F;
    new_Pup(hddm_subj(:,6)==2)=prPup_D;
    hddm_subj=[hddm_subj new_Pup];
    
     hddm_subj_probes=[hddm_subj nan(size(hddm_subj,1),8)];
    for nP=1:60
        hddm_subj_probes(hddm_subj_probes(:,4)==nP,size(hddm_subj,2)+1:size(hddm_subj_probes,2))=repmat(probe_res(nP,31:38),sum(hddm_subj_probes(:,4)==nP),1);
    end
    hddm_subj_probes(hddm_subj_probes(:,5)<-20 | hddm_subj_probes(:,5)>0,:)=[];
    table_res=[table_res ; hddm_subj_probes];
end
%% transform into tables and export
tbl_headers={'SubID','BlockN','TrialN','ProbeN','DistPr','Task','StimCat','Perf','RCode','RT'};
load(path_PsychFTlayout);
for nCh=1:63
    tbl_headers=[tbl_headers {sprintf('W_%s',layout.label{nCh})}];
end
tbl_headers=[tbl_headers ,{'Pup','pcPup','Look','State','Ori','Awa','Int','Eng','SPerf','Alert'}];
tbl_hddm=array2table(table_res,'VariableNames',tbl_headers);
tbl_hddm.State(tbl_hddm.State==4)=3;
tbl_hddm.SubID=categorical(tbl_hddm.SubID);
tbl_hddm.StimCat=categorical(tbl_hddm.StimCat);
% writetable(tbl_hddm,[root_path filesep 'hddm' filesep 'HDDM_WIM_localsleep_pup_Dec20.txt']);
tbl_hddm_F=tbl_hddm(tbl_hddm.Task==1,:);
tbl_hddm_D=tbl_hddm(tbl_hddm.Task==2,:);

tbl_hddm.Task=categorical(tbl_hddm.Task);
tbl_hddm.State=categorical(tbl_hddm.State);


%%
addpath((path_fieldtrip)); ft_defaults;

figure; %set(gcf,'position',[78         730        1162         248]);
% subplot(1,3,1)
temp_topo=mean(table2array(tbl_hddm(:,11:73)));
% temp_pV=EffectW_onPerf_pV_full(:,6)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; %caxis([-1 1]*8)
title('SW dens')
hold on;

%%
figure; %
set(gcf,'position',[78         730        1162         248]);
for nState=1:3
    subplot(2,3,nState)
    temp_topo=nanmean(table2array(tbl_hddm(tbl_hddm.State==num2str(nState),11:73)));
    % temp_pV=EffectW_onPerf_pV_full(:,6)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; caxis([0 0.2])
    title('SW dens')
    hold on;
    
    subplot(2,3,nState+3)
    temp_topo=nanmean(table2array(tbl_hddm(tbl_hddm.State==num2str(nState),11:73)))-...
        nanmean(table2array(tbl_hddm(tbl_hddm.State==num2str(1),11:73)));
    % temp_pV=EffectW_onPerf_pV_full(:,6)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; caxis([-1 1]*0.05)
    title('SW dens')
    hold on;
end
% scatter(layout.pos(find(temp_pV<0.05/63),1),...
%     layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
%     'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

%% Vigilance Ratings
figure; %
set(gcf,'position',[ 78          19        1138         959]);
for k=1:4
    subplot(4,4,k)
    temp_topo=mean(table2array(tbl_hddm(tbl_hddm.Alert==k,11:73)));
    % temp_pV=EffectW_onPerf_pV_full(:,6)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; %caxis([-1 1]*8)
    title('SW dens')
    hold on;
    caxis([0 0.2])
    
    if k>1
        subplot(4,4,k+4)
        temp_topo=mean(table2array(tbl_hddm(tbl_hddm.Alert==k,11:73)))-mean(table2array(tbl_hddm(tbl_hddm.Alert==1,11:73)));
        % temp_pV=EffectW_onPerf_pV_full(:,6)';
        simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
        colorbar; %caxis([-1 1]*8)
        title('SW dens')
        hold on;
    caxis([-1 1]*0.1)
    end
    
    subplot(4,4,k+8)
    temp_topo=mean(table2array(tbl_hddm(tbl_hddm.pcPup==k,11:73)));
    % temp_pV=EffectW_onPerf_pV_full(:,6)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; %caxis([-1 1]*8)
    title('SW dens')
    hold on;
        caxis([0 0.2])

     if k>1
        subplot(4,4,k+12)
        temp_topo=mean(table2array(tbl_hddm(tbl_hddm.pcPup==k,11:73)))-mean(table2array(tbl_hddm(tbl_hddm.pcPup==1,11:73)));
        % temp_pV=EffectW_onPerf_pV_full(:,6)';
        simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
        colorbar; %caxis([-1 1]*8)
        title('SW dens')
        hold on;
      caxis([-1 1]*0.1)
  end
end

%%
clear EffectVig_onW*
fprintf('%3.0f%%\n',round(0/(63)*100))
for nE=1:63
    fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
%     mdl_F_0= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
    mdl_F_1= fitglme(tbl_hddm,sprintf('W_%s~1+Alert+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
%     mdl_F_comp = compare(mdl_F_0,mdl_F_1);
    
        EffectVig_onW_tV(nE,:)=double(mdl_F_1.Coefficients(:,4))';
                EffectVig_onW_pV(nE,:)=double(mdl_F_1.Coefficients(:,6))';

end

figure; 
temp_topo=EffectVig_onW_tV(:,2)';
temp_pV=EffectVig_onW_pV(:,2)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('SW dens')
hold on;
scatter(layout.pos(find(temp_pV<0.05/63),1),...
    layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

%%
clear EffectVig_onW*
fprintf('%3.0f%%\n',round(0/(63)*100))
for nE=1:63
    fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
%     mdl_F_0= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
    mdl_F_1= fitglme(tbl_hddm,sprintf('W_%s~1+BlockN+TrialN+State+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
%     mdl_F_comp = compare(mdl_F_0,mdl_F_1);
    
        EffectVig_onW_tV(nE,:)=double(mdl_F_1.Coefficients(:,4))';
                EffectVig_onW_pV(nE,:)=double(mdl_F_1.Coefficients(:,6))';

end

figure; 
temp_topo=EffectVig_onW_tV(:,4)';
temp_pV=EffectVig_onW_pV(:,4)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; %caxis([-1 1]*8)
title('SW dens')
hold on;
scatter(layout.pos(find(temp_pV<0.05/63),1),...
    layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

figure; 
temp_topo=EffectVig_onW_tV(:,5)';
temp_pV=EffectVig_onW_pV(:,5)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; %caxis([-1 1]*8)
title('SW dens')
hold on;
scatter(layout.pos(find(temp_pV<0.05/63),1),...
    layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);