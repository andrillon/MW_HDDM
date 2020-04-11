%% Initilaise - clear all variables and figures
clear all
close all

%% run local file to set up paths
run ../localdef_wanderIM
addpath(genpath(lscpTools_path)) % Thomas' general toolkit
addpath(genpath(spm12_path)) % SMP12 toolbox (EEG)

data_path=[root_path filesep 'behav/']; % path of behavioural data
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
Fs_EEG=500;

check=0;
%%
hddm_res=[];
all_onsets=[];
all_ERP=[];
all_ERP2=[];
redo=1;
if redo==1
for n=1:length(files)
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
    
    if check
        %%% Load EEG
        D=spm_eeg_load([eeg_path2 filesep 'probe_infEEG_S' SubID]);
%         these_times=D.indsample(-20):D.indsample(0)-1;
        temp_data=D(1:63,:,:); % D contains the data with channels * time * trials
        temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials
    end
    
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
    load([eeg_path2 filesep 'wanderIM_twa5_noica_bigwnd_' SubID])
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs_EEG);
   fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,4)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>max_posampl| all_Waves(:,15)>max_posampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,4)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>max_posampl| all_Waves(:,15)>max_posampl,:)=[];
    
    
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
   
    %%% clean trials
    ori_clean_start_trial=clean_start_trial;
    mindistPr=nan(1,length(clean_start_trial));
    probeidx=nan(1,length(clean_start_trial));
    for nTr=1:length(clean_start_trial)
        [mindistPr(nTr),probeidx(nTr)]=findclosest(clean_start_trial(nTr)-start_probe,0);
    end
    mindistPr=mindistPr/500;
    clean_start_trial(mindistPr<-32 | mindistPr>0)=[];
    test_res(mindistPr<-32 | mindistPr>0,:)=[];
    probeidx(mindistPr<-32 | mindistPr>0)=[];
    
    temp_perf=min(test_res(:,11:12),[],2);
    temp_cat=(test_res(:,5)==test_res(:,6));
    code_resp=nan(length(temp_perf),1);
    code_resp(temp_perf==1 & temp_cat==0)=1;
    code_resp(temp_perf==1 & temp_cat==1)=0;
    code_resp(temp_perf==0 & temp_cat==0)=0;
    code_resp(temp_perf==0 & temp_cat==1)=1;
    temp_RT=(test_res(:,10)-test_res(:,8));
    
     if size(test_res,1)~=length(clean_start_trial)
        warning('... different number of trials')
    end
    % Pupil
    load([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])
    
    % Loop across trials
    fprintf('%3.0f%%\n',0)
    hddm_subj=nan(length(clean_start_trial),13+length(myElecs));
    if check
        temp_ERP=[];
        temp_ERP2=[];
    end
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
        this_probeidx=probeidx(nTr);
        this_probe=probe_res((probe_res(:,4)-1)*10+probe_res(:,1)==this_probeidx,:);
        this_MS=this_probe(32);
        if this_MS==4
            this_MS=3;
        end
                this_Vig=this_probe(38);

        %         % discard face blocks
        %         if test_res(nTr,2)==1
        %             continue;
        %         end
        
        % local sleep
        this_begTr=clean_start_trial(nTr);
        dist_toprobe=(this_begTr-start_probe(this_probeidx))/500;
        if nTr~=length(clean_start_trial) && (clean_start_trial(nTr+1)-clean_start_trial(nTr))<750
            %             this_endTr=this_begTr+500;
            this_endTr=clean_start_trial(nTr+1);
        else
            %             this_endTr=clean_start_trial(nTr+1);
            continue;
        end
        these_Waves=slow_Waves(slow_Waves(:,2)==this_probeidx,:);
        wvidx=double(these_Waves(:,5))+start_probe(this_probeidx);
        find_waves=find(wvidx>this_begTr & wvidx<this_endTr);
        this_wave=zeros(1,length(myElecs));
        this_waveF=zeros(1,length(myElecs));
        this_waveF(these_Waves(find_waves,3))=1;
%         this_waveA=nan(1,length(myElecs));
%         uniqueW=unique(these_Waves(find_waves,3));
%         for nnW=1:length(uniqueW)
%         this_waveA(uniqueW(nnW))=nanmean(these_Waves(these_Waves(find_waves,3)==uniqueW(nnW),4));
%         end
        
        if check==1
            for i=1:length(find_waves)
                if these_Waves(find_waves(i),3)==24 && these_Waves(find_waves(i),5)>-15500 && these_Waves(find_waves(i),5)<0
                    temp_eeg=temp_data(these_Waves(find_waves(i),3),these_Waves(find_waves(i),5)+32*500+(-1*D.fsample:1*D.fsample),this_probeidx);
                    temp_eeg=temp_eeg-mean(temp_eeg(1:D.fsample/2));
                    temp_ERP=[temp_ERP ; temp_eeg];
                end
            end
            
            this_begTr2=this_begTr-start_probe(this_probeidx)+32*500;
            if this_begTr2>101 && this_begTr2<15501
             temp_eeg=temp_data(47,this_begTr2+(-0.1*D.fsample:0.5*D.fsample),this_probeidx);
                    temp_eeg=temp_eeg-mean(temp_eeg(1:50));
                    temp_ERP2=[temp_ERP2 ; temp_eeg];
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
        this_pupav=nanmean(EL_data.filt_pupilSize(this_pupbegTridx:this_pupendTridx));
        hddm_subj(nTr,:)=[n this_blockidx this_trialidx this_probeidx dist_toprobe this_blockcond this_cat this_perf this_code this_rt this_waveF this_MS this_Vig this_pupav];
        
        all_onsets=[all_onsets ; [repmat([n this_blockcond this_cat],length(find_waves),1) these_Waves(find_waves,3) these_Waves(find_waves,5)+start_probe(this_probeidx)-this_begTr]];
    end
    %     hddm_res=[hddm_res ; hddm_temp(~isnan(hddm_temp(:,1)),:)];
    hddm_subj=hddm_subj(~isnan(hddm_subj(:,1)),:);
    save([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v4'],'hddm_subj')
    if check
        all_ERP=[all_ERP ; nanmean(temp_ERP(max(abs(temp_ERP),[],2)<500,:))];
        all_ERP2=[all_ERP2 ; nanmean(temp_ERP2(max(abs(temp_ERP2),[],2)<500,:))];
    end
end
end
%% Gather all individual datafiles
hddm_res=[]; nc=0;
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if ~ismember(SubID,GoodSudID)
        fprintf('skip (bad ID)\n')
        continue;
    end
    if exist([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v4.mat'])==0
         fprintf('skip (missing file)\n')
       continue;
    end
    fprintf('... %s\n',SubID)
    nc=nc+1;
    allSubID{nc}=SubID;
    load([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_v4'])
    
    hddm_subj(:,1)=str2num(SubID);
    hddm_subj(hddm_subj(:,5)<-20,:)=[];

    Pup_F=hddm_subj(hddm_subj(:,6)==1,end);
    boundaries=prctile(Pup_F,0:20:100);
    prPup_F=nan(size(Pup_F,1),1);
    for nPer=1:length(0:20:100)-2
        prPup_F(Pup_F>=boundaries(nPer) & Pup_F<boundaries(nPer+1))=nPer;
    end
        prPup_F(Pup_F>=boundaries(length(0:20:100)-1))=length(0:20:100)-1;

    Pup_D=hddm_subj(hddm_subj(:,6)==2,end);
    boundaries=prctile(Pup_D,0:20:100);
    prPup_D=nan(size(Pup_D,1),1);
    for nPer=1:length(0:20:100)-2
        prPup_D(Pup_D>=boundaries(nPer) & Pup_D<boundaries(nPer+1))=nPer;
    end
    prPup_D(Pup_D>=boundaries(length(0:20:100)-1))=length(0:20:100)-1;
    %     figure;
    %     subplot(1,2,1); plot(Pup_D); subplot(1,2,2); plot(Pup_F);
    
    new_Pup=hddm_subj(:,end);
    new_Pup(hddm_subj(:,6)==1)=prPup_F;
    new_Pup(hddm_subj(:,6)==2)=prPup_D;
    
    RT_F=hddm_subj(hddm_subj(:,6)==1,10);
    boundaries=prctile(hddm_subj(hddm_subj(:,6)==1 & hddm_subj(:,7)==0 & hddm_subj(:,10)>=0.3,10),0:10:100);
    prRT_F=nan(size(RT_F,1),1);
    for nPer=1:length(0:10:100)-2
        prRT_F(RT_F>=boundaries(nPer) & RT_F<boundaries(nPer+1))=nPer;
    end
    prRT_F(RT_F>=boundaries(length(0:10:100)-1))=length(0:10:100)-1;
        prRT_F(RT_F<0.3)=nan;

    RT_D=hddm_subj(hddm_subj(:,6)==2,10);
    boundaries=prctile(hddm_subj(hddm_subj(:,6)==2 & hddm_subj(:,7)==0 & hddm_subj(:,10)>=0.3,10),0:10:100);
    prRT_D=nan(size(RT_D,1),1);
    for nPer=1:length(0:10:100)-2
        prRT_D(RT_D>=boundaries(nPer) & RT_D<boundaries(nPer+1))=nPer;
    end
    prRT_D(RT_D>=boundaries(length(0:10:100)-1))=length(0:10:100)-1;
    prRT_D(RT_D<0.3)=nan;
    
     new_RT=hddm_subj(:,10);
    new_RT(hddm_subj(:,6)==1)=prRT_F;
    new_RT(hddm_subj(:,6)==2)=prRT_D;
    new_RT(hddm_subj(:,7)==1)=NaN;
    
    hddm_res=[hddm_res ; [hddm_subj new_Pup new_RT]];
end
%% transform into tables and export
tbl_headers={'SubID','BlockN','TrialN','ProbeN','DistProbe','Task','StimCat','Perf','RCode','RT'};
load(path_PsychFTlayout);
for nCh=1:63
     tbl_headers=[tbl_headers {sprintf('W_%s',layout.label{nCh})}];
end
tbl_headers=[tbl_headers ,{'State','Vig','Pup','pcPup','pcRT'}];
tbl_hddm=array2table(hddm_res,'VariableNames',tbl_headers);
tbl_hddm.SubID=categorical(tbl_hddm.SubID);
tbl_hddm.StimCat=categorical(tbl_hddm.StimCat);
tbl_hddm.stimulus=double(tbl_hddm.StimCat=='0');
tbl_hddm.response=double(~isnan(tbl_hddm.RT));
writetable(tbl_hddm,[root_path filesep 'hddm' filesep 'HDDM_WIM_localsleep_pup_Dec21_v5.txt']);

