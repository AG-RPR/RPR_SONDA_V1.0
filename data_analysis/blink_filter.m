% detect blinks in the time-series and remove them with
% interpolation
% new_y_resp = blink_filter(x_resp, y_resp, thresh, method)
%
% x_resp and y_resp are the raw eyetracking data
% thresh is a value (for our data, 300 recommended) that it is used to separate true
% saccades from blinking artefacts. It's in pixel/frame

function [filt_tseries, pre_fill, nanflag] = blink_filter(tseries, thresh, do_plot)

% before anything, cap the boundaries of the response in case there are
% weird artifacts
thresh_nan = length(tseries)*(1/3); % 33% of the total length
nanflag = 0;
nan_number1 = nnz(isnan(tseries));

if nan_number1 > thresh_nan
    nanflag = 1;
end

yv = diff(tseries);

zeri = find(yv==0);
tmp = diff(zeri);
zeri_index = find(tmp==1);
zeri2rmv = zeri(zeri_index);

spikes = find(abs(yv)>thresh);

to_rmv = sort([zeri2rmv spikes]);
if numel(to_rmv) == 0
    disp('No blinks found, returning time series without filtering');
    filt_tseries = tseries;
    pre_fill = tseries;
    filt_tseries = fillgaps(pre_fill,10);
    return
end
wr = 5; % window of removal (not to be confounded with window of interpolation)
for i = 1:numel(to_rmv)
    rmv_left(i,1:wr) = [to_rmv(i)-wr:to_rmv(i)-1];
    rmv_right(i,1:wr) = [to_rmv(i)+1:to_rmv(i)+wr];
end

rmv_left = reshape(rmv_left,1,numel(rmv_left));
rmv_right = reshape(rmv_right,1,numel(rmv_right));

tot_rmv = unique(sort([rmv_left rmv_right to_rmv]));
tot_rmv = tot_rmv(tot_rmv>0);
tot_rmv = tot_rmv(tot_rmv<length(tseries));

tmp2 = tseries;
tmp2(tot_rmv) = NaN;
pre_fill = tmp2;

nan_number = nnz(isnan(pre_fill));

if nan_number > thresh_nan
    nanflag = 1;
end

filt_tseries = fillgaps(pre_fill,10);

if do_plot
    subplot(2,1,1);
    plot(tseries);
    xlim([0 2400])
    title('original signal')
    subplot(2,1,2)
    plot(filt_tseries);
    xlim([0 2400])
    title('blink filtered')
    shg;
end

end
