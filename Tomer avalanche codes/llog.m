function llog(x)

if ~nargin
    tog(0)
    tog(1)
else 
    switch lower(x)
        case 'x'
            tog(1)
        case 'y'
            tog(0)
        case {'xy','yx'}
            tog(0)
            tog(1)
        otherwise
            x = str2num(x);
            if x
                set(gca,'xscale','log','yscale','log')
            else
                set(gca,'xscale','linear','yscale','linear')
            end
    end
end
shg
end

function tog(x)

if x
    sc = get(gca,'xscale');
    if strcmp(sc,'log')
        set(gca,'xscale','linear')
    else
        set(gca,'xscale','log')
    end
else
    sc = get(gca,'yscale');
    if strcmp(sc,'log')
        set(gca,'yscale','linear')
    else
        set(gca,'yscale','log')
    end
end

end