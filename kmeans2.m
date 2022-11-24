%%Reference
%Orlando Ramirez Barron (2022).k-means clustering MATLAB Central File Exchange. 
%https://www.mathworks.com/matlabcentral/fileexchange/71796-k-means-clustering
% The team has modified the open source code from the Mathwork.
clc; clear all; close all; % clear memory and command window
%%
data= readtable("adl2021.csv");
clc; clear all; close all; % clear memory and command window
%%
data= readtable("202114.csv");
%% Data of the algorithm
k=3; % Number of clusters
iteration=15; % Number of iterations
%% Generate random data points 
n=111; % number of samples or points
rx=data.Var5; % in
ry=data.Var6; % out
rz=string(data.Var4);
point=[rx,ry,rz];

%% Initial values of the centroids
cx=rx( ceil(rand(k,1)*size(rx,1)) ,:); % initial cluster centers x
cy=ry( ceil(rand(k,1)*size(ry,1)) ,:); % initial cluster centers y
%% Iterative process
for N=1:iteration
    for i=1:n
        for j=1:k
            if N==1
                distance(i,j)=sqrt(((rx(i))-(cx(j)))^2+((ry(i))-(cy(j)))^2);
            % the sample of the jth centroid
            else
                distance(i,j)=sqrt(((rx(i))-(NEWcenters_x(j)))^2+((ry(i))-(NEWcenters_y(j)))^2);
            end
        end
        % Define clusters
        [mini, CN] = min(distance(i,1:k)); % minimum distance and the
        % cluster which the sample belongs to
        Distance(N,i)=mini; Cln(N,i)=CN;
    end
    % Recompute the clusters center
    for q=1:k
        Postion_cluster=(Cln(N,:)==q); % Position of the points of the cluster
        PC(q,:)=Postion_cluster; % Points of the cluster
        NEWcenters_x(q,:)=mean(rx(Postion_cluster)); % New cluster centers in x
        NEWcenters_y(q,:)=mean(ry(Postion_cluster)); % New cluster centers in y
    end
    CLUSTER_x(N,:)=NEWcenters_x; % center of the cluster at each iteration in x
    CLUSTER_y(N,:)=NEWcenters_y; % center of the cluster at each iteration in y
    CPP(N,:,:)=PC; % points of the cluster at each iteration
end
class = zeros(1,n);
for i = 1:n
    class(i) = find(PC(:,i)==1);
end

result = table(rz, rx,ry,class');
result=sortrows(result,4);
z = result.rz;
x = result.rx;
y = result.ry;
l = result.Var4;
z=string(z);


fid = fopen('CA2021.csv', 'w');    % make csv
for i=1:n
    fprintf(fid, '%s,%d,%d,%d\n', z(i), x(i),y(i),l(i)); 
end
fclose(fid);

%% Plot of the movements of the centroids and clusters
CV= 'o*+s^v.db+c+m+k+yorobocomokoysrsbscsmsksy'; % Color Vector
for N=1:iteration+1
    figure (1)
    if N==1
        plot(rx,ry,'o','LineWidth',1.5); hold on; plot(cx,cy,'*k','LineWidth',6.5);
        hold off
    else
        for i=1:k
        plot(rx(CPP(N-1,i,:)),ry(CPP(N-1,i,:)),CV(i),'LineWidth',2); % Plot points with determined color and shape
        hold on
        end
        plot(CLUSTER_x(N-1,:),CLUSTER_y(N-1,:),'*k','LineWidth',6); hold off
         legend({'Classification','Classification','Classification','centroids'});
        title('Sydney Airport passengers first half year 2021')
        xlabel('Passengers in') 
        ylabel('Passengers out') 
    end
    grid on
    pause(0.8)
end
