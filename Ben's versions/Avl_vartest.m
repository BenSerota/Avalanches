%% Significance tests & bar plot
for i = 1:length(params)
    parameters{i} = sts(param_rows(i),:)';
    if ~LZC_flag
        data = cell(1,length(conds));
        for ii = 1:length(conds)
            data{ii} = parameters{i}(grp==ii);
        end
        if i == 1 || i == 2 % if param = sigma or alpha
            alph_sig{i} = data; % saving ala values of parameters
        end
    end
    % 1. run F test
    P = BensAnovaTest(parameters{i},alph,params{i});
    
    % 2. run paired t-tests
    if P <= alph
        [H, Pt, p_inds] = BensTtest(data,alph);
        
        % 3. correct for mult comp
        [Pt_cor, crit_p, h] = fdr_bh(Pt,alph);
        
        % 4. prepare P values for Bar graph
        Ps4bar = prepP(Pt_cor,p_inds);
        
        Param_to_bar = cellfun(@(x) mean(x),data,'UniformOutpu',false');
        Param_to_bar = cell2mat(Param_to_bar);
        save_bar = 0;
        %     E = cellfun(@(x) mean(x,2), STEs);
        %         E = cell2mat(STEs);
        E =  cell2mat(cellfun(@(x) std(x)/sqrt(length(x)), data,'UniformOutpu',false));
        
        BensSuperbar(Param_to_bar,3,Ps4bar,E,save_bar,LZC_nohb_outpath,params{i},alph)
    end
    
end

tilefigs

%% Excessory Funcs

function Ps4bar = prepP(Pt_cor,inds)
global conds
Ps4bar = zeros(length(conds));
Ps4bar(inds) = Pt_cor;
Ps4bar = Ps4bar + Ps4bar';
Ps4bar(Ps4bar==0) = 1;
end
