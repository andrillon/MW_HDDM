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

prticle_Thr_all=[99];
AmpCriterion_all=[4 9]; AmpCriterion_labels={'P2P','NP'};
thrElec_all={'byCh','allCh'};
for nThr=1:length(prticle_Thr_all)
    for nCrit=1:length(AmpCriterion_all)
        for nCrit2=1:length(AmpCriterion_labels)
            %% Parameter slow-wave detection
            prticle_Thr=prticle_Thr_all(nThr); % 80 or 90 or 95
            LimFrqW=[1 4]; % [1 4] or [4 10]
            AmpCriterionIdx=AmpCriterion_all(nCrit); % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
            thisChannel=17;
            art_ampl=200;
            max_posampl=75;
            max_Freq=7;
            fixThr=[];
            myElecs=1:63;
            Fs_EEG=500;
            
            check=0;
            %%
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
                    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>max_posampl| abs(all_Waves(:,15))>art_ampl)*100)
                    all_Waves(all_freq>max_Freq | all_Waves(:,4)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>max_posampl| abs(all_Waves(:,15))>art_ampl,:)=[];
                    
                    
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
                            if strcmp(thrElec_all{nCrit2},'byCh')
                            thr_Wave=prctile(thisE_Waves(:,AmpCriterionIdx),prticle_Thr);
                            elseif  strcmp(thrElec_all{nCrit2},'allCh')
                            thr_Wave=prctile(all_Waves(:,AmpCriterionIdx),prticle_Thr);
                            end
                        end
                        
                        slow_Waves=[slow_Waves ; (thisE_Waves(temp_p2p>thr_Wave,:))];
                        
                    end
                    
                    load([eeg_path filesep 'triggers_S' SubID])
                    
                    % Pupil
                    load([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])
                    
                    % Loop across probes
                    fprintf('%3.0f%%\n',0)
                    hddm_subj=nan(60,8+3*length(myElecs));
                    for nPr=1:length(start_probe)
                        fprintf('\b\b\b\b\b%3.0f%%\n',round(nPr/length(start_probe)*100))
                        % behav
                        this_probe=probe_res(((probe_res(:,4)-1)*10+probe_res(:,1))==nPr,:);
                        this_MS=this_probe(32);
                        if this_MS==4
                            this_MS=3;
                        end
                        this_Vig=this_probe(38);
                        this_blockidx=this_probe(4);
                        this_blockcond=this_probe(5);
                        this_probeN=this_probe(1);
                        
                        
                        % local sleep
                        this_begPr=start_probe(nPr);
                        these_Waves=slow_Waves(slow_Waves(:,2)==nPr,:);
                        wvidx=double(these_Waves(:,5));
                        wvidx2=double(these_Waves(:,7));
                        find_waves=find((wvidx/500)>-20 & (wvidx2/500)<0);
                        sel_waves=these_Waves(find_waves,:);
                        this_wave=zeros(1,length(myElecs));
                        this_waveF=zeros(1,length(myElecs));
                        this_waveA=nan(1,length(myElecs));
                        this_waveMA=nan(1,length(myElecs));
                        uniqueW=unique(these_Waves(find_waves,3));
                        for nnW=1:length(uniqueW)
                            this_waveF(uniqueW(nnW))=sum(sel_waves(:,3)==uniqueW(nnW));
                            this_waveA(uniqueW(nnW))=nanmean(sel_waves(sel_waves(:,3)==uniqueW(nnW),4));
                            this_waveMA(uniqueW(nnW))=nanmax(sel_waves(sel_waves(:,3)==uniqueW(nnW),4));
                        end
                        
                        
                        % pupil
                        pup_probeidx_all=find_trials(EL_events.Events.type,sprintf('^P%g$',this_probe(1)));
                        pup_probeidx_prev=find_trials(EL_events.Events.type(pup_probeidx_all-1),sprintf('^B%g',this_probe(4)));
                        pup_probeidx=pup_probeidx_all(pup_probeidx_prev);
                        
                        this_pupendPr=EL_events.Events.time(pup_probeidx);
                        this_pupbegPr=EL_events.Events.time(pup_probeidx)-20*EL_headers.Fs;
                        [~,this_pupbegPridx]=findclosest(EL_data.time,this_pupbegPr);
                        [~,this_pupendPridx]=findclosest(EL_data.time,this_pupendPr);
                        this_pupav=nanmedian(EL_data.filt_pupilSize(this_pupbegPridx:this_pupendPridx));
                        hddm_subj(nPr,:)=[n nPr this_blockidx this_probeN this_blockcond this_waveF this_waveA this_waveMA this_MS this_Vig this_pupav];
                        
                    end
                    %     hddm_res=[hddm_res ; hddm_temp(~isnan(hddm_temp(:,1)),:)];
                    hddm_subj=hddm_subj(~isnan(hddm_subj(:,1)),:);
                    save([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_perprobe_thr' num2str(prticle_Thr) '_' thrElec_all{nCrit2} '_' AmpCriterion_labels{nCrit} '_v5'],'hddm_subj')
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
                if exist([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_perprobe_thr' num2str(prticle_Thr) '_' thrElec_all{nCrit2} '_' AmpCriterion_labels{nCrit} '_v5.mat'])==0
                    fprintf('skip (missing file)\n')
                    continue;
                end
                fprintf('... %s\n',SubID)
                nc=nc+1;
                allSubID{nc}=SubID;
                load([root_path filesep 'behav' filesep 'HDDM_WIM' SubID '_localsleep_pup_perprobe_thr' num2str(prticle_Thr) '_' thrElec_all{nCrit2} '_' AmpCriterion_labels{nCrit} '_v5.mat'])
                
                hddm_subj(:,1)=str2num(SubID);
                
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
                
                
                hddm_res=[hddm_res ; [hddm_subj new_Pup]];
            end
            
            %% transform into tables and export
            % [n nPr this_blockidx this_probeN this_blockcond this_waveF this_waveA this_MS this_Vig this_pupav]
            tbl_headers={'SubID','ContProbeN','BlockN','ProbeN','Task'};
            load(path_PsychFTlayout);
            for nCh=1:63
                tbl_headers=[tbl_headers {sprintf('W_%s',layout.label{nCh})}];
            end
            for nCh=1:63
                tbl_headers=[tbl_headers {sprintf('A_%s',layout.label{nCh})}];
            end
            for nCh=1:63
                tbl_headers=[tbl_headers {sprintf('MA_%s',layout.label{nCh})}];
            end
            
            tbl_headers=[tbl_headers ,{'State','Vig','Pup','pcPup'}];
            tbl_hddm=array2table(hddm_res,'VariableNames',tbl_headers);
            tbl_hddm.SubID=categorical(tbl_hddm.SubID);
            writetable(tbl_hddm,[root_path filesep 'hddm' filesep 'HDDM_WIM_localsleep_amp_pup_Dec21_perprobe_thr' num2str(prticle_Thr) '_' thrElec_all{nCrit2} '_' AmpCriterion_labels{nCrit} '_v5.txt']);
        end
    end
end
