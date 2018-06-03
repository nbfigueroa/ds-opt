function xd = limitcycle_DS(x)
xd = x(2,:);
xd(2,:) = -x(1,:) + 0.9*(1 - x(1,:).^2).*x(2,:);