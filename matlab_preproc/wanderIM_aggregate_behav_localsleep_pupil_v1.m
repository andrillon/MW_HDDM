%% Init
clear all
close all

% local file, set up paths
run ../localdef_wanderIM
addpath(genpath(lscpTools_path))

data_path=[root_path filesep 'behav/'];
eeg_path=[root_path filesep 'preproc_eeg'];
pupil_path=[root_path filesep 'eyetracker'];

files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

%% Parameter slow-wave detection
prticle_Thr=80; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
thisChannel=17;

%%
hddm_res=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if exist([eeg_path filesep 'triggers_S' SubID '.mat'])==0 || exist([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])==0 || exist([eeg_path filesep 'triggers_S' SubID '.mat'])==0
        continue;
    end
     fprintf('... %s\n',SubID)
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
    RTs=[RTs; test_res(:,10)-test_res(:,8)];
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
    load([eeg_path filesep 'triggers_S' SubID])
    if size(test_res,1)~=length(clean_start_trial)
        warning('... different number of trials')
    end
    %%% first approximation: detect slow-waves on Oz (17) Pz (13) Cz (24) and
    %%% Fz (2)
    myElecs=[17 13 24 2];
    myWaves=all_Waves(ismember(all_Waves(:,3),myElecs),:);
    myWaves(double(1./((myWaves(:,7)-myWaves(:,5))/500))>4 | (myWaves(:,11)-myWaves(:,9))>150,:)=[];
    wvp2p=double(myWaves(:,11)-myWaves(:,9));
    frewv=double(1./((myWaves(:,7)-myWaves(:,5))/500));
    wvidx=double(myWaves(:,5));
    
    % Pupil
    load([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])
    
    % Loop across trials
    fprintf('%3.0f%%\n',0)
    hddm_subj=nan(length(clean_start_trial),12);
    for nTr=1:length(clean_start_trial)
        fprintf('\b\b\b\b\b%3.0f%%\n',round(nTr/length(clean_start_trial)*100))
        % behav
        this_rt=temp_RT(nTr);
        this_perf=temp_perf(nTr);
        this_code=code_resp(nTr);
        this_cat=temp_cat(nTr);
        this_blockidx=test_res(nTr,1);
        this_trialidx=test_res(nTr,4);
        
        % discard face blocks
        if test_res(nTr,2)==1
            continue;
        end
        
        % local sleep
        this_begTr=clean_start_trial(nTr);
        if nTr==length(clean_start_trial)
            this_endTr=this_begTr+500;
        else
            this_endTr=clean_start_trial(nTr+1);
        end
        find_waves=find(wvidx>this_begTr & wvidx<this_endTr);
        this_wave=zeros(1,length(myElecs));
        if ~isempty(find_waves)
            for nw=1:length(find_waves)
                if this_wave(find(ismember(myElecs,myWaves(find_waves(nw),3))))~=0
                    this_wave(find(ismember(myElecs,myWaves(find_waves(nw),3))))=max(this_wave(find(ismember(myElecs,myWaves(find_waves(nw),3)))),wvp2p(find_waves(nw)));
                else
                    this_wave(find(ismember(myElecs,myWaves(find_waves(nw),3))))=wvp2p(find_waves(nw));
                end
            end
        end
        
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
        this_pupav=nanmean(EL_data.clean_pupilSize(this_pupbegTridx:this_pupendTridx));
        hddm_subj(nTr,:)=[n this_blockidx this_trialidx this_cat this_perf this_code this_rt this_wave this_pupav];
    end
%     hddm_res=[hddm_res ; hddm_temp(~isnan(hddm_temp(:,1)),:)];
    hddm_subj=hddm_subj(~isnan(hddm_subj(:,1)),:);
    save([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup'],'hddm_subj')
end

%% Gather all individual datafiles
hddm_res=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if exist([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup.mat'])==0
        continue;
    end
       fprintf('... %s\n',SubID)
  load([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup'])
    hddm_res=[hddm_res ; hddm_subj];
end
%% transform into tables and export
tbl_headers={'SubID','BlockN','TrialN','StimCat','Perf','RCode','RT','W_Oz','W_Pz','W_Cz','W_Fz','Pup'};
tbl_hddm=array2table(hddm_res,'VariableNames',tbl_headers);
tbl_hddm.SubID=categorical(tbl_hddm.SubID);
tbl_hddm.StimCat=categorical(tbl_hddm.StimCat);
tbl_hddm.nW_Oz=tbl_hddm.W_Oz>75;
tbl_hddm.nW_Pz=tbl_hddm.W_Pz>75;
tbl_hddm.nW_Cz=tbl_hddm.W_Cz>75;
tbl_hddm.nW_Fz=tbl_hddm.W_Fz>75;

writetable(tbl_hddm,[root_path filesep 'hddm' filesep 'HDDM_WIM_localsleep_pup_May3.txt']);

mdlCorr=fitlme(tbl_hddm,'Perf~StimCat*(BlockN+TrialN+nW_Oz+nW_Pz+nW_Cz+nW_Fz+Pup)+(1|SubID)');
mdlRT=fitlme(tbl_hddm,'RT~StimCat*(BlockN+TrialN+nW_Oz+nW_Pz+nW_Cz+nW_Fz+Pup)+(1|SubID)');

%%
allRT=[];
allSW=[];
allCorr=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if exist([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup.mat'])==0
        continue;
    end
        fprintf('... %s\n',SubID)
  load([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup'])
    allCorr=[allCorr ; (hddm_subj(:,5))];
    allRT=[allRT ; nanzscore(hddm_subj(:,7))];
    allSW=[allSW ; (hddm_subj(:,9)>75)];
end
