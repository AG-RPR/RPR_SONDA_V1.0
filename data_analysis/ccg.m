function [x_axis, vel_xcor] = ccg(resp_vel, stim_vel, trim)
vel_xcor = xcorr(resp_vel,stim_vel);                                                 % cross correlation
vel_xcor = vel_xcor/norm(vel_xcor);                                                  % normalisation
vel_xcor = vel_xcor(round(end/2)-trim:round(end/2)+trim);                            % trim to +-1 sec
x_axis = -trim:trim;
end