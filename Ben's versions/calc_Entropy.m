function [] = calc_Entropy(SUBJECTS)
% this calculates Bentropy, based on -sum(PlogP)
% notice: this is based on avalanche SIZE only, though could be modified
% easily to be based on other measures

subj = length(SUBJECTS); % getting # of subjects
Size = size(SUBJECTS{1,1}.size); % getting # of thresholds and time bins 
N = nan(Size); % preallocating
temp_Bent = nan(Size); % preallocating
norm_Bent = nan(Size); % preallocating
labels = cellstr(num2str((1:subj)', 'sbj%-d')); % creating legend labels


for s = 1:subj % per subject
    for i = 1:Size(1) % different thresholds
        for j = 1:Size(2) % time bin sizes
            if i==1 && j==2 && s==2
                sprintf ('s')
            end
            % now calculating based on avalanch SIZE !!
            
            temp = full(SUBJECTS{s}.size{i,j}); % (changing sparse to full)
            temp(temp==0) = []; % eliminating zeros
            N(i,j) = round(sum(temp)); % saving N of events
            temp_Bent(i,j) = Bentropy(temp); % calculating entropy
            norm_Bent(i,j) = temp_Bent(i,j)/log2(N(i,j)); % normalizing entropy 
            %to (o,1) scale (as entropy is blocked by ln(N))
%             clear temp
        end
    end
BENT(s,:) = mean(temp_Bent,1); % average across thresholds
norm_BENT (s,:)= mean(norm_Bent,1); % average across thresholds

% clearing variables and simultaneously preallocating
N = nan(Size); % preallocating
temp_Bent = nan(Size); % preallocating
norm_Bent = nan(Size); % preallocating

end

ave_ent = mean(BENT,1); % average across subjects
ave_norm_BENT = mean(norm_BENT,1); % average across subjects

figure()
% plot([1:Size(2)], BENT); 
plot(BENT'); 
hold on
plot(ave_ent,'-og')
legend ([labels;'AVERAGE'],'location','best')
title('Entropy by bin sizes')
xlabel('time bin size')
ylabel('entropy')

% rescaling to ~(0,1)
figure()
plot(norm_BENT')
hold on
plot(ave_norm_BENT,'-og')
title('Entropy by bin sizes, rescaled ~(0,1) (normalized by logN)')
legend ([labels;'AVERAGE'],'location','best')
xlabel('time bin size')
ylabel('entropy')


function [ntrpy] = Bentropy (histog)
% just in case:
histog(histog==0) = [];       % eliminating zeros
histog = full(histog);        % sparse to double

probs = histog./sum(histog);

ntrpy = -sum(probs.*log2(probs));


