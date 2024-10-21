function offset = bestOffset(taps, grid, aimingPeriod)

x = -aimingPeriod:0.001:aimingPeriod;

for i = 1:length(x)
diff(i) = sum(abs(grid-taps+x(i)));
end



[minvalue, idx] = min(diff);



offset = x(idx);
end