function [MEEG, frame] = onesuj(ntri,ntime,effects);


MEEG = randn(ntri,ntime);

frame = [];
for i = 1:size(MEEG,1)
    if i < ntri/2
        frame(i).cond = 'blue';
    else
        frame(i).cond = 'red';
    end
    frame(i).dist =randn;
end

effect = strcmp({frame.cond},'red') .* effects(1) + [frame.dist] .* effects(2) + [frame.dist] .*  strcmp({frame.cond},'red') .* effects(3);

MEEG = bsxfun(@plus,MEEG,effect');
