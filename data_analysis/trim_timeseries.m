function [stim, resp]=trim_timeseries(stim, resp)
% trim the first and last samples of each time series to avoid possible artifacts in recordings 
    stim = stim(60:end-60);
    resp = resp(60:end-60);
    
    [cuthere, who2cut] = min([length(stim) length(resp)]);
    switch who2cut
        case 1
            resp = resp(1:cuthere);
        case 2
            stim = stim(1:cuthere);
    end
end