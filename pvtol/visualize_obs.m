function visualize_obs(obs)

ellipse(obs(:,3),obs(:,3),zeros(size(obs,1),1),obs(:,1),obs(:,2),'k',[],1);
axis equal