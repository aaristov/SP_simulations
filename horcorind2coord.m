% in - index
% out - um using FOV 50 um
function out = horcorind2coord(in) %(xp/FOV+1/2)*512+32
out=((in-33)/512-1/2)*50;
end