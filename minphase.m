function [minphpilots] = minphase(K, Lp, L)
movratio = 0.1;
pilot = randn(K*L, Lp);
for lk = 1 : K*L
    pilot(lk, :) = pilot(lk, :) / norm(pilot(lk, :));
end
cor = abs(pilot * pilot');
[corr, rowin] = max(cor-eye(K*L));
[maxcor, colin] = max(corr);
rowin = rowin(1, colin);
maxcorpre = maxcor * 2;
count = 0;
counts = 0;
while counts < 1e1
    maxcorpre = maxcor;
    pilot1 = pilot(rowin, :);
    pilot2 = pilot(colin, :);
    pilotd = (pilot1 - pilot2) * sign(pilot1 * pilot2');
    pilot11 = pilot1 + pilotd * movratio;
    pilot22 = pilot2 - pilotd * movratio;
    pilot(rowin, :) = pilot11 / norm(pilot11);
    pilot(colin, :) = pilot22 / norm(pilot22);
    cor = abs(pilot * pilot');
    [corr, rowin] = max(cor-eye(K*L));
    [maxcor, colin] = max(corr);
    rowin = rowin(1, colin);
    count = count + 1;
    if maxcor-maxcorpre > 0
        counts = counts + 1;
    end
end
minphpilots = [pilot(1:K,:),pilot(K+1:2*K,:),pilot(2*K+1:3*K,:)];
