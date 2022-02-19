function visualize_obs(obs,color)
if nargin<2
    color = 'k';
end
ellipse(obs(:,3),obs(:,3),zeros(size(obs,1),1),obs(:,1),obs(:,2),color,[],1);
axis equal