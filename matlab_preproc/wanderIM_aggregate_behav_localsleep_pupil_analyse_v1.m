%% Init
clear all
close all

% local file, set up paths
run ../localdef_wanderIM
addpath(genpath(lscpTools_path))

data_path=[root_path filesep 'behav/'];
eeg_path=[root_path filesep 'preproc_eeg'];
pupil_path=[root_path filesep 'eyetracker'];
eeg_path2=[root_path filesep 'preproc_ica'];

files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

load([pwd filesep '..' filesep '..' filesep 'WanderIM' filesep 'paper' filesep 'paper_SubID'])

%% Parameter slow-wave detection
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
thisChannel=17;
art_ampl=200;
max_posampl=75;
max_Freq=7;
fixThr=[];
myElecs=1:63;
%%
hddm_res=[];
for n=7:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if exist([eeg_path filesep 'triggers_S' SubID '.mat'])==0 || exist([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])==0 || exist([eeg_path filesep 'triggers_S' SubID '.mat'])==0
        continue;
    end
    fprintf('... %s\n',SubID)
    if ~ismember(SubID,GoodSudID)
        continue;
    end
    %     CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
    % SART
    %  1: num block
    %  2: block cond (1: Faces / 2: Squares)
    %  3: image set
    %  4: num trial
    %  5: seq trial
    %  6: target
    %  7: resp
    %  8: stim onset
    %  9: stim pre
    % 10: resp onset
    % 11: nogo
    % 12: go
    
    %%% Gathering behavioural data
    %     RTs=[RTs; test_res(:,10)-test_res(:,8)];
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(n,nbt)=nanmean(tp_gos);
        corr_nogo(n,nbt)=nanmean(tp_nogos);
    end
    temp_perf=min(test_res(:,11:12),[],2);
    temp_cat=(test_res(:,5)==test_res(:,6));
    code_resp=nan(length(temp_perf),1);
    code_resp(temp_perf==1 & temp_cat==0)=1;
    code_resp(temp_perf==1 & temp_cat==1)=0;
    code_resp(temp_perf==0 & temp_cat==0)=0;
    code_resp(temp_perf==0 & temp_cat==1)=1;
    temp_RT=(test_res(:,10)-test_res(:,8));
    
    % Local sleep
    load([eeg_path filesep 'MWCTR_cont_twa2_' SubID])
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./500);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,2)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,2)>art_ampl | all_Waves(:,11)>max_posampl,:)=[];
    
    slow_Waves=[];
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        if ~isempty(fixThr)
            thr_Wave=fixThr;
        else
            thr_Wave=prctile(all_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        
        slow_Waves=[slow_Waves ; (thisE_Waves(temp_p2p>thr_Wave,:))];
        
    end
    
    load([eeg_path filesep 'triggers_S' SubID])
    if size(test_res,1)~=length(clean_start_trial)
        warning('... different number of trials')
    end
    
    
    % Pupil
    load([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])
    
    % Loop across trials
    fprintf('%3.0f%%\n',0)
    hddm_subj=nan(length(clean_start_trial),9+length(myElecs));
    for nTr=1:length(clean_start_trial)
        fprintf('\b\b\b\b\b%3.0f%%\n',round(nTr/length(clean_start_trial)*100))
        % behav
        this_rt=temp_RT(nTr);
        this_perf=temp_perf(nTr);
        this_code=code_resp(nTr);
        this_cat=temp_cat(nTr);
        this_blockidx=test_res(nTr,1);
        this_blockcond=test_res(nTr,2);
        this_trialidx=test_res(nTr,4);
        
        %         % discard face blocks
        %         if test_res(nTr,2)==1
        %             continue;
        %         end
        
        % local sleep
        this_begTr=clean_start_trial(nTr);
        if nTr==length(clean_start_trial)
            this_endTr=this_begTr+500;
        else
            this_endTr=clean_start_trial(nTr+1);
        end
        wvidx=double(slow_Waves(:,5));
        find_waves=find(wvidx>this_begTr & wvidx<this_endTr);
        this_waveF=zeros(1,length(myElecs));
        this_waveF(slow_Waves(find_waves,3))=1;
        
        
        % pupil
        pup_trialidx=find_trials(EL_events.Events.type,sprintf('^B%g_T%g$',this_blockidx,this_trialidx));
        pup_trialidxnext=find_trials(EL_events.Events.type,sprintf('^B%g_T%g$',this_blockidx,this_trialidx+1));
        this_pupbegTr=EL_events.Events.time(pup_trialidx);
        if isempty(pup_trialidxnext)
            this_pupendTr=this_pupbegTr+EL_headers.Fs;
        else
            this_pupendTr=EL_events.Events.time(pup_trialidxnext);
        end
        [~,this_pupbegTridx]=findclosest(EL_data.time,this_pupbegTr);
        [~,this_pupendTridx]=findclosest(EL_data.time,this_pupbegTr);
        this_pupav=nanmean(EL_data.filt_pupilSize(this_pupbegTridx:this_pupendTridx));
        hddm_subj(nTr,:)=[n this_blockidx this_trialidx this_blockcond this_cat this_perf this_code this_rt this_waveF this_pupav];
    end
    %     hddm_res=[hddm_res ; hddm_temp(~isnan(hddm_temp(:,1)),:)];
    hddm_subj=hddm_subj(~isnan(hddm_subj(:,1)),:);
    save([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v2'],'hddm_subj')
end

%% Gather all individual datafiles
hddm_res=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if ~ismember(SubID,GoodSudID)
        continue;
    end
    if exist([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v2.mat'])==0
        continue;
    end
    fprintf('... %s\n',SubID)
    load([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v2'])
    Pup_F=hddm_subj(hddm_subj(:,4)==1,end);
    boundaries=prctile(Pup_F,0:20:100);
    prPup_F=nan(size(Pup_F,1),1);
    for nPer=1:5
        prPup_F(Pup_F>=boundaries(nPer) & Pup_F<boundaries(nPer+1))=nPer;
    end
    
    Pup_D=hddm_subj(hddm_subj(:,4)==2,end);
    boundaries=prctile(Pup_D,0:20:100);
    prPup_D=nan(size(Pup_D,1),1);
    for nPer=1:5
        prPup_D(Pup_D>=boundaries(nPer) & Pup_D<boundaries(nPer+1))=nPer;
    end
    
    %     figure;
    %     subplot(1,2,1); plot(Pup_D); subplot(1,2,2); plot(Pup_F);
    
    new_Pup=hddm_subj(:,end);
    new_Pup(hddm_subj(:,4)==1)=prPup_F;
    new_Pup(hddm_subj(:,4)==2)=prPup_D;
    hddm_res=[hddm_res ; [hddm_subj new_Pup]];
end
%% transform into tables and export
tbl_headers={'SubID','BlockN','TrialN','Task','StimCat','Perf','RCode','RT'};
load(path_PsychFTlayout);
for nCh=1:63
    tbl_headers=[tbl_headers {sprintf('W_%s',layout.label{nCh})}];
end
tbl_headers=[tbl_headers ,{'Pup','pcPup'}];
tbl_hddm=array2table(hddm_res,'VariableNames',tbl_headers);
tbl_hddm.SubID=categorical(tbl_hddm.SubID);
tbl_hddm.StimCat=categorical(tbl_hddm.StimCat);
writetable(tbl_hddm,[root_path filesep 'hddm' filesep 'HDDM_WIM_localsleep_pup_Dec20.txt']);
tbl_hddm_F=tbl_hddm(tbl_hddm.Task==1,:);
tbl_hddm_D=tbl_hddm(tbl_hddm.Task==2,:);

tbl_hddm.Task=categorical(tbl_hddm.Task);

%% LMEs on Waves
fprintf('%3.0f%%\n',0)
EffectTask_onW=[];
for nE=1:63
    fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
    %     mdl_F_0= fitlme(tbl_hddm,sprintf('W_%s~1+BlockN+TrialN+(1|SubID)',layout.label{nE}));
    %     mdl_F_1= fitlme(tbl_hddm,sprintf('W_%s~1+Task*(BlockN+TrialN)+(1|SubID)',layout.label{nE}));
    %     mdl_F_comp = compare(mdl_F_0,mdl_F_1);
    
    mdl_F_0= fitglme(tbl_hddm,sprintf('W_%s~1+BlockN+TrialN+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm,sprintf('W_%s~1+Task*(BlockN+TrialN)+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
        EffectTask_onW(nE,:)=[double(mdl_F_1.Coefficients(4,[4 6])) double(mdl_F_1.Coefficients(5,[4 6])) double(mdl_F_1.Coefficients(6,[4 6])) NaN double(mdl_F_1.Coefficients(2,[4 6])) double(mdl_F_1.Coefficients(3,[4 6]))];
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
        EffectTask_onW(nE,:)=[double(mdl_F_1.Coefficients(4,[4 6])) double(mdl_F_1.Coefficients(5,[4 6])) double(mdl_F_1.Coefficients(6,[4 6])) double(mdl_F_comp.pValue(2)) double(mdl_F_1.Coefficients(2,[4 6])) double(mdl_F_1.Coefficients(3,[4 6]))];
    end
    
    EffectTask_onW(nE,:)=[double(mdl_F_1.Coefficients(4,[4 6])) double(mdl_F_1.Coefficients(5,[4 6])) double(mdl_F_1.Coefficients(6,[4 6])) double(mdl_F_comp.pValue(2)) double(mdl_F_1.Coefficients(2,[4 6])) double(mdl_F_1.Coefficients(3,[4 6]))];
    
end

%%
addpath((path_fieldtrip)); ft_defaults;
figure; set(gcf,'position',[78         730        1162         248]);
subplot(1,2,1)
temp_topo=EffectTask_onW(:,8)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Block')
hold on;
scatter(layout.pos(find(EffectTask_onW(:,9)<fdr(EffectTask_onW(:,9),0.05)),1),...
    layout.pos(find(EffectTask_onW(:,9)<fdr(EffectTask_onW(:,9),0.05)),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

subplot(1,2,2)
temp_topo=EffectTask_onW(:,10)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Trial')
hold on;
scatter(layout.pos(find(EffectTask_onW(:,11)<fdr(EffectTask_onW(:,11),0.05)),1),...
    layout.pos(find(EffectTask_onW(:,11)<fdr(EffectTask_onW(:,11),0.05)),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);


figure; set(gcf,'position',[78         730        1162         248]);
subplot(1,3,1)
temp_topo=EffectTask_onW(:,1)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Task on SW')
hold on;
scatter(layout.pos(find(EffectTask_onW(:,2)<0.05/63),1),...
    layout.pos(find(EffectTask_onW(:,2)<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

subplot(1,3,2)
temp_topo=EffectTask_onW(:,3)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Task*Block')
hold on;
scatter(layout.pos(find(EffectTask_onW(:,4)<0.05/63),1),...
    layout.pos(find(EffectTask_onW(:,4)<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

subplot(1,3,3)
temp_topo=EffectTask_onW(:,5)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Task*Trial')
hold on;
scatter(layout.pos(find(EffectTask_onW(:,6)<0.05/63),1),...
    layout.pos(find(EffectTask_onW(:,6)<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);
%% LMEs by electrodes
%
% fprintf('%3.0f%%\n',0)
% for nE=1:63
%     fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
%     mdl_F_0= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
%     mdl_F_1= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
%     mdl_F_comp = compare(mdl_F_0,mdl_F_1);
%
% %     EffectW_onGoRT_F(nE,:)=[double(mdl_F_1.Coefficients(4,2)) double(mdl_F_1.Coefficients(4,4)) double(mdl_F_1.Coefficients(4,6)) double(mdl_F_comp.pValue(2))];
%
%
%     mdl_F_0= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+(1|SubID)'));
%     mdl_F_1= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
%     mdl_F_comp = compare(mdl_F_0,mdl_F_1);
%
% %     EffectW_onGo_F(nE,:)=[double(mdl_F_1.Coefficients(4,2)) double(mdl_F_1.Coefficients(4,4)) double(mdl_F_1.Coefficients(4,6)) double(mdl_F_comp.pValue(2))];
%
%
%     mdl_F_0= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='1',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
%     mdl_F_1= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='1',:),sprintf('RT~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
%     mdl_F_comp = compare(mdl_F_0,mdl_F_1);
%
% %     EffectW_onNoGo_F(nE,:)=[double(mdl_F_1.Coefficients(4,2)) double(mdl_F_1.Coefficients(4,4)) double(mdl_F_1.Coefficients(4,6)) double(mdl_F_comp.pValue(2))];
%
%     %%%%%% DIGITS
%     mdl_D_0= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
%     mdl_D_1= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
%     mdl_D_comp = compare(mdl_D_0,mdl_D_1);
%
% %     EffectW_onGoRT_D(nE,:)=[double(mdl_D_1.Coefficients(4,2)) double(mdl_D_1.Coefficients(4,4)) double(mdl_D_1.Coefficients(4,6)) double(mdl_D_comp.pValue(2))];
%
%
%     mdl_D_0= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+(1|SubID)'));
%     mdl_D_1= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
%     mdl_D_comp = compare(mdl_D_0,mdl_D_1);
%
% %     EffectW_onGo_D(nE,:)=[double(mdl_D_1.Coefficients(4,2)) double(mdl_D_1.Coefficients(4,4)) double(mdl_D_1.Coefficients(4,6)) double(mdl_D_comp.pValue(2))];
%
%
%     mdl_D_0= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='1',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
%     mdl_D_1= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='1',:),sprintf('RT~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
%     mdl_D_comp = compare(mdl_D_0,mdl_D_1);
%
% %     EffectW_onNoGo_D(nE,:)=[double(mdl_D_1.Coefficients(4,2)) double(mdl_D_1.Coefficients(4,4)) double(mdl_D_1.Coefficients(4,6)) double(mdl_D_comp.pValue(2))];
% end

%% LMEs by electrodes

fprintf('%3.0f%%\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
    mdl_F_0= fitlme(tbl_hddm(tbl_hddm.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+Task+(1|SubID)'));
    mdl_F_1= fitlme(tbl_hddm(tbl_hddm.StimCat=='0',:),sprintf('RT~(1+BlockN+TrialN+Task)*W_%s+(1|SubID)',layout.label{nE}));
    mdl_F_comp = compare(mdl_F_0,mdl_F_1);
    
    EffectW_onGoRT_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onGoRT_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onGoRT_chi(nE)=mdl_F_comp.pValue(1);
    
    
    mdl_F_0= fitglme(tbl_hddm(tbl_hddm.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+Task+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm(tbl_hddm.StimCat=='0',:),sprintf('Perf~(1+BlockN+TrialN+Task)*W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onGo_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onGo_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onGo_chi(nE)=mdl_F_comp.pValue(1);
    
    mdl_F_0= fitglme(tbl_hddm(tbl_hddm.StimCat=='1',:),sprintf('Perf~1+BlockN+TrialN+Task+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm(tbl_hddm.StimCat=='1',:),sprintf('Perf~(1+BlockN+TrialN+Task)*W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onNoGo_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onNoGo_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onNoGo_chi(nE)=mdl_F_comp.pValue(1);
end

%%
clear EffectW_*_F* EffectW_*_D*
fprintf('%3.0f%%\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
    mdl_F_0= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
    mdl_F_1= fitlme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
    mdl_F_comp = compare(mdl_F_0,mdl_F_1);
    
    EffectW_onGoRT_F_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onGoRT_F_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onGoRT_F_chi(nE)=mdl_F_comp.pValue(1);
    
    
    mdl_F_0= fitglme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm_F(tbl_hddm_F.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onGo_F_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onGo_F_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onGo_F_chi(nE)=mdl_F_comp.pValue(1);
    
    mdl_F_0= fitglme(tbl_hddm_F(tbl_hddm_F.StimCat=='1',:),sprintf('Perf~1+BlockN+TrialN+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm_F(tbl_hddm_F.StimCat=='1',:),sprintf('Perf~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onNoGo_F_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onNoGo_F_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onNoGo_F_chi(nE)=mdl_F_comp.pValue(1);
    
    mdl_F_0= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+(1|SubID)'));
    mdl_F_1= fitlme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('RT~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}));
    mdl_F_comp = compare(mdl_F_0,mdl_F_1);
    
    EffectW_onGoRT_D_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onGoRT_D_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onGoRT_D_chi(nE)=mdl_F_comp.pValue(1);
    
    
    mdl_F_0= fitglme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm_D(tbl_hddm_D.StimCat=='0',:),sprintf('Perf~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onGo_D_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onGo_D_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onGo_D_chi(nE)=mdl_F_comp.pValue(1);
    
    mdl_F_0= fitglme(tbl_hddm_D(tbl_hddm_D.StimCat=='1',:),sprintf('Perf~1+BlockN+TrialN+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm_D(tbl_hddm_D.StimCat=='1',:),sprintf('Perf~1+BlockN+TrialN+W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onNoGo_D_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onNoGo_D_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onNoGo_D_chi(nE)=mdl_F_comp.pValue(1);
end

%%
figure;
lineLabels={'GoRT','GoPerf','NogoPerf'};
for nline=1:3
    if nline==1
        tvals_F=EffectW_onGoRT_F_tV;
        pvals_F=EffectW_onGoRT_F_pV;
        
        tvals_D=EffectW_onGoRT_D_tV;
        pvals_D=EffectW_onGoRT_D_pV;
        
        tvals=EffectW_onGoRT_tV;
        pvals=EffectW_onGoRT_pV;
    elseif nline==2
        tvals_F=EffectW_onGo_F_tV;
        pvals_F=EffectW_onGo_F_pV;
        
        tvals_D=EffectW_onGo_D_tV;
        pvals_D=EffectW_onGo_D_pV;
        
        tvals=EffectW_onGo_tV;
        pvals=EffectW_onGo_pV;
    elseif nline==3
        tvals_F=EffectW_onNoGo_F_tV;
        pvals_F=EffectW_onNoGo_F_pV;
        
        tvals_D=EffectW_onNoGo_D_tV;
        pvals_D=EffectW_onNoGo_D_pV;
        
        tvals=EffectW_onNoGo_tV;
        pvals=EffectW_onNoGo_pV;
    end
    subplot(3,3,3*(nline-1)+1);
    temp_topo=tvals_F(:,4)';
    temp_pV=pvals_F(:,4)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; caxis([-1 1]*8)
    title([lineLabels{nline} ' - Faces'])
    hold on;
    scatter(layout.pos(find(temp_pV<0.05/63),1),...
        layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
        'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);
    
    subplot(3,3,3*(nline-1)+2);
    temp_topo=tvals_D(:,4)';
    temp_pV=pvals_D(:,4)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; caxis([-1 1]*8)
    title([lineLabels{nline} ' - Digits'])
    hold on;
    scatter(layout.pos(find(temp_pV<0.05/63),1),...
        layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
        'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);
    
    subplot(3,3,3*(nline-1)+3);
    temp_topo=tvals(:,8)';
    temp_pV=pvals(:,8)';
    simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
    colorbar; caxis([-1 1]*8)
    title([lineLabels{nline} '*Task'])
    hold on;
    scatter(layout.pos(find(temp_pV<0.05/63),1),...
        layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
        'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);
end

%%

fprintf('%3.0f%%\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
    %     mdl_F_0= fitglme(tbl_hddm,sprintf('Perf~1+StimCat*Task*BlockN*TrialN+(1|SubID)'),'Distribution','Binomial','Link','logit');
    %     mdl_F_1= fitglme(tbl_hddm,sprintf('Perf~(1+StimCat*Task*BlockN*TrialN)*W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    mdl_F_0= fitglme(tbl_hddm,sprintf('Perf~1+StimCat*Task+(1|SubID)'),'Distribution','Binomial','Link','logit');
    mdl_F_1= fitglme(tbl_hddm,sprintf('Perf~(1+StimCat*Task)*W_%s+(1|SubID)',layout.label{nE}),'Distribution','Binomial','Link','logit');
    %     mdl_F_2= fitlme(tbl_hddm,sprintf('Perf~(1+StimCat*Task)*W_%s+(1|SubID)',layout.label{nE}));
    if mdl_F_0.LogLikelihood>mdl_F_1.LogLikelihood
        mdl_F_comp =compare(mdl_F_1,mdl_F_0);
    else
        mdl_F_comp =compare(mdl_F_0,mdl_F_1);
    end
    EffectW_onPerf_tV_full(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
    EffectW_onPerf_pV_full(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
    EffectW_onPerf_chi_full(nE)=mdl_F_comp.pValue(1);
end

%%
figure; set(gcf,'position',[78         730        1162         248]);
subplot(1,3,1)
temp_topo=EffectW_onPerf_tV_full(:,6)';
temp_pV=EffectW_onPerf_pV_full(:,6)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Task on SW')
hold on;
scatter(layout.pos(find(temp_pV<0.05/63),1),...
    layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

subplot(1,3,2)
temp_topo=EffectW_onPerf_tV_full(:,16)';
temp_pV=EffectW_onPerf_pV_full(:,16)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Task on SW')
hold on;
scatter(layout.pos(find(temp_pV<0.05/63),1),...
    layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

subplot(1,3,3)
temp_topo=EffectW_onPerf_tV_full(:,26)';
temp_pV=EffectW_onPerf_pV_full(:,26)';
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'off',[],0,1);
colorbar; caxis([-1 1]*8)
title('Task on SW')
hold on;
scatter(layout.pos(find(temp_pV<0.05/63),1),...
    layout.pos(find(temp_pV<0.05/63),2),'Marker','o','MarkerFaceColor',[1 1 1]*0.7,...
    'MarkerEdgeColor',[1 1 1]*0.7,'SizeData',24);

%% HITS and FAs
uniqueSub=unique(tbl_hddm.SubID);
tbl_Hits=[];
tbl_Misses=[];
tbl_FAs=[];
for nS=1:length(uniqueSub)
    for nB=1:6
        temp_tbl=tbl_hddm(tbl_hddm.SubID==uniqueSub(nS) & tbl_hddm.BlockN==nB,:);
        Hits=mean(temp_tbl.RCode(temp_tbl.StimCat=='0')==1);
        Misses=mean(temp_tbl.RCode(temp_tbl.StimCat=='0')==0);
        FAs=mean(temp_tbl.RCode(temp_tbl.StimCat=='1')==1);
        
        RT_Hits=mean(temp_tbl.RT(temp_tbl.StimCat=='0' & temp_tbl.RCode==1));
        RT_FAs=mean(temp_tbl.RT(temp_tbl.StimCat=='1' & temp_tbl.RCode==1));
        
        physio=mean(table2array(temp_tbl(:,9:73)),1);
        physio_Hits=mean(table2array(temp_tbl(temp_tbl.StimCat=='0' & temp_tbl.RCode==1,9:73)),1);
        physio_Misses=mean(table2array(temp_tbl(temp_tbl.StimCat=='0' & temp_tbl.RCode==0,9:73)),1);
        physio_FAs=mean(table2array(temp_tbl(temp_tbl.StimCat=='1' & temp_tbl.RCode==1,9:73)),1);
        temp_mat1=[nS nB double(unique(temp_tbl.Task)) Hits RT_Hits physio_Hits];
        temp_mat2=[nS nB double(unique(temp_tbl.Task)) Misses nan(size(Misses)) physio_Misses];
        temp_mat3=[nS nB double(unique(temp_tbl.Task)) FAs RT_FAs physio_FAs];
        
        tbl_Hits=[tbl_Hits ; temp_mat1];
        tbl_Misses=[tbl_Misses ; temp_mat2];
        tbl_FAs=[tbl_FAs ; temp_mat3];
    end
end
tbl_headers={'SubID','BlockN','Task','Perf','RT'};
for nCh=1:63
    tbl_headers=[tbl_headers {sprintf('W_%s',layout.label{nCh})}];
end
tbl_headers=[tbl_headers ,{'Pup','pcPup'}];
tbl_Hits=array2table(tbl_Hits,'VariableNames',tbl_headers);
tbl_Misses=array2table(tbl_Misses,'VariableNames',tbl_headers);
tbl_FAs=array2table(tbl_FAs,'VariableNames',tbl_headers);

% %%
% clear EffectW_onHits*
% fprintf('%3.0f%%\n',0)
% for nE=1:63
%     fprintf('\b\b\b\b\b%3.0f%%\n',round(nE/(63)*100))
%     mdl_F_0= fitlme(tbl_FAs,sprintf('Perf~1+Task+(1|SubID)'));
%     mdl_F_1= fitlme(tbl_FAs,sprintf('Perf~1+Task*W_%s+(1|SubID)',layout.label{nE}));
%     mdl_F_comp =compare(mdl_F_0,mdl_F_1);
%     
%     EffectW_onHits_tV(nE,:)=[double(mdl_F_1.Coefficients(:,4))'];
%     EffectW_onHits_pV(nE,:)=[double(mdl_F_1.Coefficients(:,6))'];
%     EffectW_onHits_chi(nE)=mdl_F_comp.pValue(1);
% end


%%
uniqueSub=unique(tbl_hddm.SubID);
for nS=1:length(uniqueSub)
    fprintf('%g/%g\n',nS,length(uniqueSub))
    for nT=1:2
        temp_tbl=tbl_hddm(tbl_hddm.SubID==uniqueSub(nS) & tbl_hddm.Task==num2str(nT),:);
        for nE=1:63
            [bb,dev,stats] = glmfit(temp_tbl.Perf(temp_tbl.StimCat=='0'),table2array(temp_tbl(temp_tbl.StimCat=='0',match_str(temp_tbl.Properties.VariableNames,['W_' layout.label{nE}]))));
            beta_GoResp(nS,nT,nE)=stats.beta(2);
            beta_pV_GoResp(nS,nT,nE)=stats.p(2);
            
            [bb,dev,stats] = glmfit(temp_tbl.Perf(temp_tbl.StimCat=='1'),table2array(temp_tbl(temp_tbl.StimCat=='1',match_str(temp_tbl.Properties.VariableNames,['W_' layout.label{nE}]))));
            beta_NoGoResp(nS,nT,nE)=stats.beta(2);
            beta_pV_NoGoResp(nS,nT,nE)=stats.p(2);
        end
    end
end

for nE=1:63
    for nT=1:2
        [bb,dev,stats] = glmfit(tbl_hddm.Perf(tbl_hddm.StimCat=='0' & tbl_hddm.Task==num2str(nT)),table2array(tbl_hddm(tbl_hddm.StimCat=='0' & tbl_hddm.Task==num2str(nT),match_str(tbl_hddm.Properties.VariableNames,['W_' layout.label{nE}]))));
        beta_GoResp2(nT,nE)=stats.beta(2);
        beta_pV_GoResp2(nT,nE)=stats.p(2);
        
        [bb,dev,stats] = glmfit(tbl_hddm.Perf(tbl_hddm.StimCat=='1' & tbl_hddm.Task==num2str(nT)),table2array(tbl_hddm(tbl_hddm.StimCat=='1' & tbl_hddm.Task==num2str(nT),match_str(tbl_hddm.Properties.VariableNames,['W_' layout.label{nE}]))));
        beta_NoGoResp2(nT,nE)=stats.beta(2);
        beta_pV_NoGoResp2(nT,nE)=stats.p(2);
    end
    
    [bb,dev,stats] = glmfit([tbl_hddm.Perf(tbl_hddm.StimCat=='0') tbl_hddm.Task(tbl_hddm.StimCat=='0')=='2'],table2array(tbl_hddm(tbl_hddm.StimCat=='0',match_str(tbl_hddm.Properties.VariableNames,['W_' layout.label{nE}]))));
    beta_GoResp3(nE)=stats.beta(2);
    beta_pV_GoResp3(nE)=stats.p(2);
    
    [bb,dev,stats] = glmfit([tbl_hddm.Perf(tbl_hddm.StimCat=='1') tbl_hddm.Task(tbl_hddm.StimCat=='1')=='2'],table2array(tbl_hddm(tbl_hddm.StimCat=='1',match_str(tbl_hddm.Properties.VariableNames,['W_' layout.label{nE}]))));
    beta_NoGoResp3(nE)=stats.beta(2);
    beta_pV_NoGoResp3(nE)=stats.p(2);
end
%%
figure; format_fig;
[~,hb(1)]=simpleTplot(1:63,squeeze(beta_GoResp(:,1,:)),0,'r',[3 0.05 0],'-',0.5,1,[],[],2); hold on;
[~,hb(2)]=simpleTplot(1:63,squeeze(beta_GoResp(:,2,:)),0,'b',[3 0.05 0],'-',0.5,1,[],[],2); hold on;

figure; format_fig;
plot(1:63,squeeze(beta_GoResp2(1,:)),'Color','r','LineWidth',2); hold on;
plot(1:63,squeeze(beta_GoResp2(2,:)),'Color','b','LineWidth',2); hold on;

figure; format_fig;
plot(1:63,squeeze(beta_GoResp3),'Color','r','LineWidth',2); hold on;


%%
figure; format_fig;
[~,hb(1)]=simpleTplot(1:63,squeeze(beta_NoGoResp(:,1,:)),0,'r',[3 0.05 0],'-',0.5,1,[],[],2); hold on;
[~,hb(2)]=simpleTplot(1:63,squeeze(beta_NoGoResp(:,2,:)),0,'b',[3 0.05 0],'-',0.5,1,[],[],2); hold on;

figure; format_fig;
plot(1:63,squeeze(beta_NoGoResp2(1,:)),'Color','r','LineWidth',2); hold on;
plot(1:63,squeeze(beta_NoGoResp2(2,:)),'Color','b','LineWidth',2); hold on;

%%
figure; format_fig;
plot(1:63,zscore(beta_GoResp3),'Color','r','LineWidth',2); hold on;
plot(1:63,zscore(beta_NoGoResp3),'Color','b','LineWidth',2); hold on;
