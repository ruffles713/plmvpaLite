function [subj results] = JR_run_mvpa_v3(mvpa_workspace)

tic

%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flags.save_workspace = 1; % 1 = yes please, 0 = no thanks

% unless a previously created workspace is specified, load in the preprocessed data
if ~exist('mvpa_workspace')

    subj_id = 's106';
    exp_name = 'PAST';
    roi_file = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/merged_AAL_ROIs.nii';
    roi_name = 'merged_AAL_ROIs'
    num_runs = 10;
    num_TP_per_run = 203;

    % load user-created filename and onsets lists into workspace
    load('raw_filenames_s8mm_wa.mat')
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);
    vol_info = spm_vol(raw_filenames{1}); %get functional data resolution info for spm .img writing
    [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);

    if flags.save_workspace == 1
        save_cmd = ['save ' subj_id '_' roi_name '_8mm_smoothing_' mvpa_datetime '.mat'];
        eval(save_cmd);
    end
else
    eval(['load ' mvpa_workspace])
end

% Set flags (% 1 = yes please, 0 = no thanks)
flags.equate_number_of_old_new_trials_per_subjective_bin = 0;
flags.equate_number_of_trials_in_cond_1_and_2 = 1;
flags.plot_mean_timecourses = 0;
flags.plot_ROC_curve = 0;
flags.display_performance_breakdown = 1;
flags.generate_importance_maps = 1;


num_results_iter = 2; % number of times to run the cross validation process

anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)

% specify which conditions to use for classification (must correspond to the names of conditions specified below)
condnames =  {'Subjective_old','Subjective_new'};

% classifier parameters
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract info about conditions from onsets file
num_conds = size(onsets,2);
all_regs = zeros(num_conds,num_runs*num_TP_per_run); % initialize regs matrix as conditions x timepoints

for cond = 1: num_conds-1 %(exclude last condition ("no_response" trials)
    for trial = 1: length(onsets{cond})
        time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
        all_regs(cond,time_idx) = 1;
    end
end


% SPECIAL FIX FOR s104 and s105 because reg #6 (NEW_recollect) is a fake
% placeholder trial
if ~isempty(find(strcmp(subj_id,{'s104','s105'})))
    all_regs(6,:)=0;
end

% Artificially balance the number of trials given each subjective memory response
if flags.equate_number_of_old_new_trials_per_subjective_bin == 1
    for j = 1:5
        OLD_trials = find(all_regs(j,:));
        NEW_trials = find(all_regs(j+5,:));

        num_OLD = length(OLD_trials);
        num_NEW = length(NEW_trials);

        if num_OLD > num_NEW
            rand_array = rand(1,num_OLD);
            [sorted inds]= sort(rand_array);
            trials_to_cut = OLD_trials(inds(1:num_OLD-num_NEW));
            all_regs(j,trials_to_cut)=0;
        elseif num_OLD < num_NEW
            rand_array = rand(1,num_NEW);
            [sorted inds]= sort(rand_array);
            trials_to_cut = NEW_trials(inds(1:num_NEW-num_OLD));
            all_regs(j+5,trials_to_cut)=0;
        end
    end
end


% define conditions of interest
% specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'

Objective_new = sum(all_regs(1:5,:));
Objective_new = sum(all_regs(6:10,:));

Subjective_old = sum(all_regs([1 2 3 6 7 8],:));
Subjective_new = sum(all_regs([4 5 9 10],:));

Subjective_old_HC_only = sum(all_regs([1 2 6 7],:));
Subjective_new_HC_only = sum(all_regs([5 10],:));

Hits = sum(all_regs([1 2 3],:));
Misses = sum(all_regs([4 5],:));
CRs = sum(all_regs([9 10],:));
FAs = sum(all_regs(6:8,:));

R_hits = all_regs(1,:);
HC_hits = all_regs(2,:);
LC_hits = all_regs(3,:);
LC_misses = all_regs(4,:);
HC_misses = all_regs(5,:);
R_FAs = all_regs(6,:);
HC_FAs = all_regs(7,:);
LC_FAs = all_regs(8,:);
LC_CRs = all_regs(9,:);
HC_CRs = all_regs(10,:);

R_and_HC_hits = sum(all_regs([1 2],:));

%no_resp = all_regs(11,:); %excluded from analysis

%assign conditions to train/test classifier on
regs_of_interest = [];
eval(['regs_of_interest(1,:) = ' condnames{1} ';'])
eval(['regs_of_interest(2,:) = ' condnames{2} ';'])


if flags.equate_number_of_trials_in_cond_1_and_2 == 1

    cond1_trials = find(regs_of_interest(1,:));
    cond2_trials = find(regs_of_interest(2,:));
    num_cond1 = length(cond1_trials);
    num_cond2 = length(cond2_trials);

    if num_cond1 > num_cond2
        rand_array = rand(1,num_cond1);
        [sorted inds]= sort(rand_array);
        trials_to_cut = cond1_trials(inds(1:num_cond1-num_cond2));
        regs_of_interest(1,trials_to_cut) = 0;
        display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{1}]);
    elseif num_cond1 < num_cond2
        rand_array = rand(1,num_cond2);
        [sorted inds]= sort(rand_array);
        trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
        regs_of_interest(2,trials_to_cut) = 0;
        display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{2}]);
    else
        display('Trial numbers are already balanced');
    end
end

display([num2str(count(regs_of_interest(1,:)==1)) ' trials in condition ' condnames{1}])
display([num2str(count(regs_of_interest(2,:)==1)) ' trials in condition ' condnames{2}])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select TRs of interest (to correspond with peak post-stim BOLD response)

all_trials = sum(all_regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)

TR0 = subj.patterns{end}.mat(:,find(all_trials)-1); % -1 TR (-2-0 sec)
TR1 = subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
TR2 = subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
TR3 = subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
TR4 = subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
TR5 = subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)

mean_data = (TR3+TR4)/2; %take uniform average of peak signal TRs
%mean_data = (TR3*.25+TR4+TR5*.25)/1.5;

% BASELINE CORRECT MEAN DATA (optional)
%mean_baseline = TR1;
%mean_data = mean_data - mean_baseline;

clear TR0 TR1 TR2 TR3 TR4 TR5; %clean up matlab workspace to save memory


if flags.plot_mean_timecourses == 1 % plot mean timecourses for whole ROI (useful check of data quality & overall activity pattern)

    c1_ind = find(regs_of_interest(1,:));
    for t = 1: length(c1_ind)
        c1_ts(t,:) = mean(subj.patterns{end}.mat(:,c1_ind(t):c1_ind(t)+6),1);
    end

    c2_ind = find(regs_of_interest(2,:));
    for t = 1: length(c2_ind)
        c2_ts(t,:) = mean(subj.patterns{end}.mat(:,c2_ind(t):c2_ind(t)+6),1);
    end

    figure;
    plot(mean(c1_ts,1),'b');
    hold on;
    plot(mean(c2_ts,1),'r');
end


% condense regs by removing zeros
trial_counter = 1;
for i = 1: size(all_regs,2)
    if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
        condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
        condensed_regs_all(:,trial_counter) = all_regs(:,i);
        condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
        trial_counter = trial_counter + 1;
    end
end

% initialize regressors object
subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% add new condensed activation pattern
subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',mean_data);

zhist = sprintf('Pattern ''%s'' created by JR custom code','spiral_d_hp_z_condensed');
subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);

% clean up workspace
subj = remove_mat(subj,'pattern','spiral_d_hp_z');
clear mean_data;

% update run vector to condensed format
subj.selectors{1}.mat = condensed_runs;
subj.selectors{1}.matsize = size(condensed_runs);

% "activate" only those trials of interest (from regs_of_interest) before
% creating cross-validation indices

nonactive_trials = find(sum(condensed_regs_of_interest)==0);
active_trials = find(sum(condensed_regs_of_interest));
actives_selector = ones(1,size(condensed_regs_all,2)); % intialize vector of all ones
actives_selector(nonactive_trials) = 0; % remove all non-"regs_of_interst" trials (set to zero)
subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
subj = set_mat(subj,'selector','conditions_of_interest',actives_selector);
subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');


% run feature selection ANOVA (if desired)
if anova_p_thresh ~= 1
    subj = feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',anova_p_thresh);
    classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
else
    classifier_mask = subj.masks{1}.name; % use original mask
end

%%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for x = 1: num_results_iter

    [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',classifier_mask,class_args);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % do some important RAM clean-up and data archiving
    for y = 1:num_runs
        if flags.generate_importance_maps == 1
            results_IW{x}.iterations(y).scratchpad.net.IW{1} = single(results{x}.iterations(y).scratchpad.net.IW{1}); % save weights to pass to JR_interpret_weights
        end
        results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
    end

    if flags.display_performance_breakdown == 1
        % analyze the results in more detail
        correct_vector = [];
        desireds_vector = [];
        guesses_vector = [];
        acts_diff_vector = [];
        for a = 1:num_runs
            correct_vector = horzcat(correct_vector,results{x}.iterations(a).perfmet.corrects);
            desireds_vector = horzcat(desireds_vector,results{x}.iterations(a).perfmet.desireds);
            guesses_vector = horzcat(guesses_vector,results{x}.iterations(a).perfmet.guesses);
            acts_diff_vector = horzcat(acts_diff_vector, results{x}.iterations(a).acts(1,:)-results{x}.iterations(a).acts(2,:));
        end

        overall_accuracy = mean(correct_vector)
        overall_hit_rate = mean(correct_vector(desireds_vector==1))
        overall_fa_rate = 1-mean(correct_vector(desireds_vector==2))
        overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate)


        % classification_accuracy_by_resp (for full classification only)
        for b = 1:10
            classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_regs_all(b,active_trials))));
            mean_acts_diffs_by_resp(b) = mean(acts_diff_vector(find(condensed_regs_all(b,active_trials))));
            number_of_trials_per_bin(b) = length(find(condensed_regs_all(b,active_trials)));
        end

        display(classification_accuracy_by_resp);
        display(mean_acts_diffs_by_resp);
        display(number_of_trials_per_bin);

        [sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
        correct_sorted = correct_vector(ind);
        desireds_sorted = desireds_vector(ind);

        % sort by absolute value of classifier "confidence"
        [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
        abs_correct_sorted = correct_vector(abs_ind);
        abs_desireds_sorted = desireds_vector(abs_ind);

        % print results of top 50,100, 150, etc.
        num_trials = length(abs_correct_sorted);
        if num_trials>50

            for j= 1:floor(num_trials/50)
                sorted_acc(j)=mean(abs_correct_sorted(1:50*j));
            end
            sorted_acc(j+1)=mean(abs_correct_sorted(1:end));

            display(sorted_acc)
        end
    end

    if flags.plot_ROC_curve == 1

        % create continuous ROC function
        for i = 1:length(sorted_diffs);
            hit_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
            fa_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
        end

        % created binned ROC
        counter = 1;
        for i = 1:8
            hit(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:counter+45]))) / length(find(desireds_sorted == 1));
            fa(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:counter+45]))) / length(find(desireds_sorted == 2));
            counter = counter+46;
        end
        figure
        plot(fa_rate,hit_rate,'.-')
        hold on
        plot([0 1],[0 1],'r')
        xlabel('P(Old|New)')
        ylabel('P(Old|Old)')
    end
end


if flags.generate_importance_maps == 1;

    subj = JR_interpret_weights(subj, results);

    immap1 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
    immap2 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2

    if anova_p_thresh == 1 % NO ANOVA VERSION

        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        immap1(voxel_inds)=subj.patterns{end}.mat(:,1); %store immap values at appropriate voxel indices
        immap2(voxel_inds)=subj.patterns{end}.mat(:,2); %store immap values at appropriate voxel indices

        immap1 = immap1*1000;
        immap2 = immap2*1000;

        vol_info.fname = [condnames{1} '_no_anova_' mvpa_datetime '.img'];
        spm_write_vol(vol_info,immap1);
        vol_info.fname = [condnames{2} '_no_anova_' mvpa_datetime '.img'];
        spm_write_vol(vol_info,immap2);

    else % ANOVA VERSION

        for j = 1:num_runs
            temp1 = zeros(vol_info.dim); %initialize appropriately sized matrix
            temp2 = zeros(vol_info.dim); %initialize appropriately sized matrix
            voxel_inds{j} = find(subj.masks{end-num_runs+j}.mat); %get mask voxel indices
            temp1(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,1); %store immap values at appropriate voxel indices
            temp2(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,2); %store immap values at appropriate voxel indices
            immap1 = immap1+temp1; %add values cumulatively across iterations
            immap2 = immap2+temp2;
        end

        %sum across masks to get composite mask (where value of each voxel =
        %number of runs for which that voxel was included)
        composite_mask = zeros(vol_info.dim);
        for i = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
            composite_mask = composite_mask+subj.masks{i}.mat;
        end
        voxels_to_exclude = find(composite_mask<6);  % exclude voxels that exist for less than 6 of the ANOVA masks
        immap1(voxels_to_exclude)=0;
        immap2(voxels_to_exclude)=0;

        immap1_avg = immap1./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
        immap2_avg = immap2./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling


        vol_info.fname = [condnames{1} '_p' num2str(anova_p_thresh) '_' mvpa_datetime '.img'];
        spm_write_vol(vol_info,immap1_avg);
        vol_info.fname = [condnames{2} '_p' num2str(anova_p_thresh) '_' mvpa_datetime '.img'];
        spm_write_vol(vol_info,immap2_avg);
    end

end

time2finish = toc/60;
display(['Finished in ' num2str(time2finish) ' minutes']);