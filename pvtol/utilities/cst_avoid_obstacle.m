%% ----- formulate the state constraints for collision avoidance 
function [c,ceq,cgrad,ceqgrad] = cst_avoid_obstacle(x,obstacles,tubeXZ_size,nu)
n_obs = size(obstacles,1);
N = size(x,2);
nx = size(x,1); 
ndec = nx+nu+1;         % time(1), 6 states, 2 inputs. 
c = zeros(n_obs,N);     % c<=0 is the path inequality constraint

for i=1:n_obs
    obs_x = obstacles(i,1);
    obs_z = obstacles(i,2);
    obs_r = obstacles(i,3);
    c(i,:) = (obs_r+tubeXZ_size).^2-((x(1,:)-obs_x).^2+(x(2,:)-obs_z).^2);
end
if nargout == 4 % analytic gradients    
    % Matlab convention
%     cgrad =  zeros(ndec*N,n_obs*N); 
%     for j=1:N
%         cgrad((j-1)*ndec+(1:2),(i-1)*N+j) = [-2*(x(1,j)-obs_x); -2*(x(2,j)-obs_z)];
%     end  
    
    % optTraj convention
    cgrad =  zeros(n_obs,ndec,N); 
    for i=1:n_obs
        obs_x = obstacles(i,1);
        obs_z = obstacles(i,2);
        for j=1:N
    %         cgrad((j-1)*ndec+(1:2),(i-1)*N+j) = [-2*(x(1,j)-obs_x); -2*(x(2,j)-obs_z)];
            cgrad(i,2:3,j) = [-2*(x(1,j)-obs_x) -2*(x(2,j)-obs_z)]; % note that the first column is t
        end
    end
    ceq = []    ;% ceq=0 is the path equality constraint
    ceqgrad = [];
    % gradc = sparse(gradc);
end
end
