clear;clc;close all;tic
A=[-100,-100;-100,100;100,-100;100,100];
theata = 0:(2*pi/720) :2*pi;
R_theata = (abs(cos(3*theata/2)) +abs(sin(3*theata/2))).^(1/2);
% CN_di(:,1) = 16.8* R_theata .* cos(theata);%%x
% CN_di(:,2) = 16.8* R_theata .* sin(theata);%%y
% 
% CW_di(:,1) = (16.8 * R_theata + 0.005) .* cos(theata);
% CW_di(:,2) = (16.8 * R_theata + 0.005) .* sin(theata);

% temp = CN_di(:,1);
% CN_di(:,1) = -CN_di(:,2); x = -y
% CN_di(:,2) = temp;
% 
% temp = CW_di(:,1);
% CW_di(:,1) = -CW_di(:,2);
% CW_di(:,2) = temp;
CN_di(:,1) = -16.8* R_theata .* sin(theata);%%x
CN_di(:,2) = 16.8* R_theata .* cos(theata);%%y

CW_di(:,1) = -(16.8 * R_theata + 0.005) .* sin(theata);
CW_di(:,2) = (16.8 * R_theata + 0.005) .* cos(theata);

pfix = [A;CN_di;CW_di];

figure(2)
hold on;
box on;
plot(CN_di(:,1),CN_di(:,2), CW_di(:,1),CW_di(:,2), 'b', 'LineWidth', 1.5);
plot(0,0,'.k', 'MarkerSize', 5);
axis equal;
ylim([-100, 100]);
xlim([-100, 100]);
set(gca, 'XTick',(-100:20:100)) % 设置x坐标轴的刻度 
set(gca, 'YTick',(-100:20:100))  % 设置y坐标轴的刻度
set(gca,'xtick',[],'ytick',[])
set(gca, 'FontSize',15, 'FontWeight', 'bold','LineWidth', 1.5) % 设置坐标轴字体是 8

%
%%
figure(1)
hold on;
%plot(CN_di(:,1), CN_di(:,2),'r',CW_di(:,1), CW_di(:,2),'b');

fd = @(p) drectangle(p,-100,100,-100,100);
fh = @(p) (0.26 + 0.3*abs( sqrt( p(:,1).^2 + p(:,2).^2 ) - 16.8*(abs(cos(3/2*atan(-p(:,1) ./ p(:,2)))) + abs(sin(3/2*atan(-p(:,1)./p(:,2))))).^(1/2) ));
[p,t] =  distmesh2d(fd,fh,0.19,[-100,-100;100,100],pfix);
% patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] );
% patch('Faces', CN_di, 'faces', t, 'EdgeColor', 'blue', 'FaceColor', 'none', 'LineWidth', '1.2' );


%癌样细胞的直角坐标方程
 %d = sqrt( p(:,1).^2 + p(:,2).^2 ) - 16.8*(abs(cos(3/2*atan(p(:,2) ./ p(:,1)))) + abs(sin(3/2*atan(p(:,2)./p(:,1))))).^(1/2);
% save(p.mat);save(cm_n.mat);save(cm_w.mat);



%%%%%%%%%%%%%%%%%%%%%%%%%%%  step2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;tic;
load("matlab1_di.mat");
[VX,VY] = voronoi(p(:,1),p(:,2));%生成Voronoi图,返回 Voronoi边的二维顶点。
h = plot(VX,VY,'-b',p(:,1),p(:,2),'.r');%绘图

boundary=[-100,100];
xlim(boundary)
ylim(boundary)

%%
%%%%%%%%%%%%%%%%%显示凸包序号%%%%%%%%%
%% MATLAB官方文档 https://ww2.mathworks.cn/help/matlab/math/voronoi-diagrams.html?searchHighlight=voronoi&s_tid=srchtitle_voronoi_2
nump = size(p,1);
plabels = arrayfun(@(n) {sprintf('p%d', n)}, (1:nump)');
hold on
Hpl = text(p(:,1), p(:,2), plabels, 'color', 'r', ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
%%
%%%%%%%%%%%%%%%%% 计算voronoi单元  %%%%%%%%%%%%%%

dt = delaunayTriangulation(p); %%基于P中的点创建二维Delaunay三角剖分     %以delaunay为基础计算voronoi参数
[V,R] = voronoiDiagram(dt);    %%返回Delaunay三角剖分中点的Voronoi顶点V和Voronoi区域r,
                               %%r 中的每个区域表示围绕某个三角剖分顶点的点，它们比三角剖分中的其他顶点更靠近该顶点

[row,col]=find(abs(V)>100);    %%find 查找非零元素的索引和值  ----绝对值>100就超出了边界
 

for i=1:size(R,1)              %%循环长度  size(R,1):返回第一列的长度                               
   if sum(R{i}==row)~=0        %% R{i}给出与点位i相关联的Voronoi顶点的索引判断有没有超出边界的点，如果和不等于0就表明有超过边界的点；    
      del=find(sum(R{i}==row));
        R{i}(del)=[];          %% 然后把这些点置为空数组。
  end
end

a1=size(V,1);                  %% a1=数组V的长度，即维诺图中单元的顶点个数
%%
hori=find(round(abs(p(:,1)),-4)==100 & p(:,2)<100 & p(:,2)>=-100);  
%%p的第一列取绝对值四舍五入，找到-50<=p第二列<50,p第一列=50 的点 （x=+-50,-50<=y<50）
%%hori表示边界框的左右边界

for i=1:size(hori,1)
    m=max(V(R{hori(i)},2));   %%取V的第二列的最大值 
    V(size(V,1)+1,:)=[p(hori(i),1),m];
end  %%这段作用？？？？
%%
vert=find(round(abs(p(:,2)),-4)==100 & p(:,1)<100 & p(:,1)>=-100);
%%vert表示边界框的上下边界(-50<=x<50,y=+-50)

for i=1:size(vert,1)
    m=max(V(R{vert(i)},1));
    V(size(V,1)+1,:)=[m,p(vert(i),2)];
end
%%VC单元排序
a2=size(V,1);
Adara(:,1)=(a1:a2)';Adara(:,2:3)=V(a1:a2,:);
Bound=find(round(abs(p(:,2)),-4)==100 | round(abs(p(:,1)),-4)==100);

for i=1:size(Bound,1)
    Adara(:,4)=sqrt((Adara(:,2)-p(Bound(i),1)).^2+(Adara(:,3)-p(Bound(i),2)).^2);
    Adara=sortrows(Adara,4);
    R{Bound(i)}=[R{Bound(i)},Adara(1,1),Adara(2,1)];
end

%%
%%%%%%%%%%%%%显示voronoi单元顶点序号%%%%%%%%%%%%
numv = size(V,1);
vlabels = arrayfun(@(n) {sprintf('V%d', n)}, (2:numv)');
hold on
Hpl = text(V(2:end,1), V(2:end,2), vlabels, ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
hold off
%%数据处理  
%%依照方框大小而定，单元间的距离不会超过方框范围，但剖分会出现inf值，需取成+-50；
for i=1:length(V)
    if V(i,1) > 100 && V(i,1) ~= Inf
       V(i,1) = 100;
    end
    if V(i,1) < -100 && V(i,1) ~= Inf
       V(i,1) = -100;
    end
    if V(i,2) > 100 && V(i,2) ~= Inf
       V(i,2) = 100;
    end
    if V(i,2) < -100 && V(i,2) ~= Inf
       V(i,2) = -100;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%  step3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=roundn(p,-10); 
V=roundn(V,-10);

%%绘制维诺图
[VX,VY] = voronoi(p(:,1),p(:,2));
h = plot(VX,VY,'-b',p(:,1),p(:,2),'.r');
xlim([-100,100]);
ylim([-100,100]);

%%
%%%%%%%%%%%邻接矩阵%%%%%%%%%%%
Adj=cell(length(p),1);
Boun=zeros(length(p));
for i=1:size(p,1)  %x 遍历
    k=1;           %%索引编号
    for j=1:size(p,1) %y 遍历
        if size(intersect(R{i},R{j}),2)==2   %%intersect 设置两个数组R{1},R{2}的交集 找出i=j的点，（对角线） 返回第二列的长度
            A=intersect(R{i},R{j});
            Adj{i,1}(k)=j;
            Boun(i,j)=sqrt((V(A(1),1)-V(A(2),1))^2+(V(A(1),2)-V(A(2),2))^2);%%计算两点间的距离
            if Boun(i,j)==0    %%如果邻接点间的距离为0即这是同一个点  两个单元的边连接一点点就忽略掉
                Adj{i,1}(k)=[];%%则把这个点的邻接矩阵中所表示的值置为0
            elseif Boun(i,j) == Inf
                Boun(i,j) = 100;
            elseif Boun(i,j) == -Inf
                Boun(i,j) = -100;
            else
                k=k+1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  step4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=roundn(p,-10);
CN_di=roundn(CN_di,-10);
CW_di=roundn(CW_di,-10);
%%Adj.mat存储邻接矩阵，Boun.mat存储点间距离
Th=10;          %%模型深度 取胞体中心至边缘最远处的距离 即R(pi/6)微米
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始条件%%%%%%%%%%%%%%%%%%%%%%%%
%
% sigma_ele  1.2  电介质电导率
% sigma_mem  9.5*10^(-9)     膜电导率
% epse epsi  电介质介电常数
% epsm       膜介电常数
%%
Mem_d = 0.005;  %%膜厚5nm
sigma_ele=1.2;
sigma_mem=9.5*1e-9;
epslong_ele=7.08*1e-10; %电介质介电常数
epslong_mem=4.43*1e-11; %细胞质介电常数

%%
%%cm 在膜上的P点的索引列向量
CN_index_di = zeros(length(CN_di),1);
for i = 1:length(p)
   for j=1:length(CN_di)
    if p(i,:) == CN_di(j,:)
        CN_index_di(i,1) = i;
    end
   end
end
CN_index_di(CN_index_di == 0) =[];

CW_index_di = zeros(length(CW_di),1);
for i = 1:length(p)
   for j=1:length(CW_di)
    if p(i,:) == CW_di(j,:)
        CW_index_di(i,1) = i;
    end
   end
end
CW_index_di(CW_index_di == 0) =[];
%%计算膜上点的索引 --- 便于下面给区域分配介电参数和电导率
%%%CM
Cell_index_di =zeros(length(CN_index_di),2);
for i=1:size(CN_index_di,1)
    for j=1:size(CW_index_di,1)
        D=sqrt((p(CN_index_di(i,1),1)-p(CW_index_di(j,1),1))^2+(p(CN_index_di(i,1),2)-p(CW_index_di(j,1),2))^2);%% 求D膜内外两点间的距离即等于膜的厚度5nm
        if roundn(D,-4) == Mem_d
            Cell_index_di(i,1)=CN_index_di(i);
            Cell_index_di(i,2)=CW_index_di(j);%Mem n*2的数组 第一列为内膜点索引，第二列为外膜点索引 //(i,j)是一对在膜上的节点
            break
        end
    end
end

%%
%%%%%%%%%%%%%%%%%分配导电系数和介电常数%%%%%%%%%%%%%%%%%%%%%
sig=zeros(size(p,1),1);
eps=zeros(size(p,1),1);

for i=1:size(p,1)
    if (size(find(Cell_index_di == i),1) == 1)
        sig(i,1) = sigma_mem;
        eps(i,1) = epslong_mem;
    else
        sig(i,1) = sigma_ele;
        eps(i,1) = epslong_ele;
    end
end

%%
%%%%%%%%%%%%%为VC单元分配电阻值电容值%%%%%%%%%%%%%%%%%%%%%
%%Y---->导纳
Res_VC=zeros(size(Boun));Cap_VC=zeros(size(Boun));
for i = 1:size(Adj,1)
    for j = 1:size(Adj{i},2)
        sigma_av = (sig(i) + sig(Adj{i}(j)))/2;
        epslong_av = (eps(i) + eps(Adj{i}(j)))/2;
        l_ij=sqrt((p(i,1)-p(Adj{i}(j),1))^2+(p(i,2)-p(Adj{i}(j),2))^2);
        Cap_VC(i,Adj{i}(j)) =(epslong_av * (Boun(i,Adj{i}(j)) * 1e-6) * (Th * 1e-6))/(l_ij * 1e-6);
        Res_VC(i,Adj{i}(j)) = (sigma_av * (Boun(i,Adj{i}(j)) * 1e-6) * (Th * 1e-6))/(l_ij * 1e-6) ;

    end
end

%%
%%%%%%%%%%%%%节点导纳矩阵%%%%%%%%%%%%%%%%%%%%%%%

k2=find(roundn(abs(p(:,2)),-4)~=100); %%K1是上下边界 K2是内部
Res=zeros(size(k2,1));
Cap=zeros(size(k2,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%《推导正负号》%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(k2,1)
    for j=1:size(k2,1)
        if i==j    %%i=j 主对角线上，求自导纳
            Res(i,j)=sum(Res_VC(k2(i),:));
            Cap(i,j)=-sum(Cap_VC(k2(i),:)); %%+-号的意义？节点电压法推导
        else     %%求互导纳 加-号
            Res(i,j)=-Res_VC(k2(i),k2(j));
            Cap(i,j)=Cap_VC(k2(i),k2(j));
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%电源导纳矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=find(roundn(p(:,2),-4)==100);k2=find(roundn(abs(p(:,2)),-4)~=100);
Rus=zeros(size(k2,1),1);
for i=1:size(k2,1)           %判断Res1中第k2行第k1列元素是否为0,
    for j=1:size(k1,1)
        if Res_VC(k2(i),k1(j))~=0
            Rus(i)=Rus(i) - Res_VC(k2(i),k1(j));  %%排除膜上临近节点电流的影响
        end
    end
end
%%

disp('完成！step_4');
toc