function viz_point(time, evo_time, X, Y)
% Point over time
figure
plot(time,squeeze(evo_time(Y,X,:)));
end

