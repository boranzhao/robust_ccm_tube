function obs = generate_obs(num_obs,xrange,yrange,rrange)

obs = zeros(num_obs,3);

% xs = (xrange(2)-xrange(1)).*rand(num_obs,1) + xrange(1);
% ys = (yrange(2)-yrange(1)).*rand(num_obs,1) + yrange(1);
% rs = (rrange(2)-rrange(1)).*rand(num_obs,1)*10 + rrange(1);
% 
% for i=1:num_obs
%     obs(i,:) = [xs(i) ys(i) rs(i)];
% end
% 


allX = 1000;
allY = 1000;
allR = 1000;
while length(allX)< num_obs+1
    thisX = (xrange(2)-xrange(1)).*rand(1,1) + xrange(1);
    thisY = (yrange(2)-yrange(1)).*rand(1,1) + yrange(1);
    thisR = (rrange(2)-rrange(1)).*rand(1,1)+ rrange(1);
    distances = sqrt((allX - thisX).^2 + (allY - thisY) .^ 2);
    if min(distances) > 2*rrange(2)+0.1
      % Nothing else is close, so add this one
      allX(end+1,1) = thisX;
      allY(end+1,1) = thisY;
      allR(end+1,1) = thisR;
    end
end
obs = [allX allY allR];
obs = obs(2:end,:);