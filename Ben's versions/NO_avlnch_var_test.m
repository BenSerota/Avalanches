%avlnch var test
new_cond_marker = [avprms.cond];
new_cond_marker(bads) = [];

%% creating final structures
a = struct();
for i = 1:length(params)
    for ii = 1:length(conds)
        a = setfield(a,params{i},conds{ii},...
            sts(param_rows(i),new_cond_marker==ii)');
    end
end

flds = {'Parameter','P_Anov','P_corrected', 'H'};
A = cell2struct(cell(length(params),length(flds)),flds,2);

% 1. run F test, per parameter
% Big_param_cell = cell(1,length(params)); % preallocating for large cross-run comparison
for i = 1:length(params)
    A(i).Parameter = params{i};
    
    % turning data into cell for accecibility
    temp = cell(1,4);
    for ii = 1:4
        temp{ii} = a.(params{i}).(conds{ii});
    end
    A(i).P_Anov = BensAnovaTest(temp,alph,params{i});
    
    if A(i).P_Anov <= alph
        % 2. run paired t-tests
        [H, Pt, p_inds] = BensTtest(temp,alph);
        
        % 3. correct for mult comp
        [A(i).P_corrected, crit_p, A(i).H] = fdr_bh(Pt,alph);
        
        % % 4. prepare P values for Bar graph
        Ps4bar = prepP(A(i).P_corrected,p_inds);
        
        Param_to_bar = cellfun(@(x) mean(x),temp,'UniformOutpu',false);
        Param_to_bar = cell2mat(Param_to_bar);
        save_bar = 0;
        %     E = cellfun(@(x) mean(x,2), STEs);
        E =  cell2mat(cellfun(@(x) std(x)/sqrt(length(x)), temp,'UniformOutpu',false));
        
        BensSuperbar(Param_to_bar,3,Ps4bar,E,save_bar,avlnch_rslts,params{i},alph)
    end
%     Big_param_cell{i}(run_count) = Param_to_bar;
    tilefigs
end

function Ps4bar = prepP(Pt_cor,inds)
global conds
Ps4bar = zeros(length(conds));
Ps4bar(inds) = Pt_cor;
Ps4bar = Ps4bar + Ps4bar';
Ps4bar(Ps4bar==0) = 1;
end