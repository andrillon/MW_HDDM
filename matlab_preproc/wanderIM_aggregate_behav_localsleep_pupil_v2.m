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
