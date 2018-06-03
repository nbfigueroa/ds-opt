function xd = globally_stable_DS(x)
xd = - x(1,:);
xd(2,:) = - x(1,:) .* cos(x(1,:)) - x(2,:);