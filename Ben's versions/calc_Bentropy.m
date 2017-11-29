function [] = calc_Bentropy(SUBJECTS)
Size = size(SUBJECTS{1,1}.size);
sub = size(SUBJECTS,1);
n = nan(Size); % preallocating
temp_BENT = nan(Size); % preallocating
labels = cellstr(num2str((1:sub)', 'sbj#%-d'));


for s = 1:sub
    for i = 1:size(SUBJECTS{s}.size,1) % 5 different thresholds
        for j = 1:size(SUBJECTS{s}.size,2) % time bin sizes
            temp = full(SUBJECTS{s}.size{i,j}); % (changing sparse to full)
            temp(temp==0) = []; % eliminating zeros
            n(i,j) = sum(temp); % saving N of events
            temp_BENT(i,j) = Bentropy(temp); % calculating entropy
            norm_Bent(i,j) = temp_BENT(i,j)/log2(n(i,j)); % normalizing entropy 
            %to (o,1) scale (as entropy is blocked by ln(N))
        end
    end
BENT(s,:) = mean(temp_BENT,1); % average across thresholds
norm_BENT (s,:)= mean(norm_Bent,1); % average across thresholds
end

ave_ent = mean(BENT,1); % for average line
ave_norm_BENT = mean(norm_BENT,1); % for average line

figure()
plot([1:j], BENT); 
hold on
plot(ave_ent,'-og')
legend ([labels;'AVERAGE'],'location','best')
title('Eentropy by bin sizes')
xlabel('time bin size')
ylabel('entropy')


% plotting into a sigmoid
siggi = 1./(1+exp(-BENT));
siggi_ave_ent = mean(siggi,1); % calculating mean (after sig transformaiton, as it's non linear)
figure()
plot(siggi')
hold on
plot(siggi_ave_ent,'-og')
title('Eentropy by bin sizes, rescaled ~(0,1) (normalized by Sigmoid function)')
legend ([labels;'AVERAGE'],'location','best')
xlabel('time bin size')
ylabel('entropy')

% rescaling to ~(0,1)
figure()
plot(norm_BENT')
hold on
plot(ave_norm_BENT,'-og')
title('Eentropy by bin sizes, rescaled ~(0,1) (normalized by logN)')
legend ([labels;'AVERAGE'],'location','best')
xlabel('time bin size')
ylabel('entropy')

