clc;
clear all
%=========== 读取每天的密集度nc文件 =======================
%读取第一天的影像
path = 'D:\wzx\amsr2\'; 
dstr = '20190315';
ncfile = [path, 'asi-AMSR2-n3125-',dstr,'-v5.4.nc'];
if ~exist(ncfile,'file')
    hdffile = [path, 'asi-AMSR2-n3125-',dstr,'-v5.4.hdf'];
    sitc1 = hdfread(hdffile, '/ASI Ice Concentration');
elseif exist(ncfile,'file')
    sitc = ncread(ncfile, 'z');  
    sitc1=sitc';   % nc文件和hdf文件的图是横纵坐标相反的，需要转置，才能和经纬度矩阵匹配
end
% figure; 
% fig = imshow(sitc1,[0,100]); colormap jet; colorbar % 显示图像
sitc1 = sitc1(1450:1850,1600:1900);
% sitc1 = sitc1(1500:1750,1650:1800);
% figure; imshow(sitc1,[-5 100]);  colormap jet; colorbar;

% ========================================================
% 读取第二天的影像
% clc;
% clear all;
path = 'D:\wzx\amsr2\';
dstr2 = '20190323';
hdffile2 = [path, 'asi-AMSR2-n3125-',dstr2,'-v5.4.hdf'];
sitc2 = hdfread(hdffile2, '/ASI Ice Concentration');
% figure; 
% fig = imshow(sitc2,[0,100]); colormap jet; colorbar % 显示图像
sitc2 = sitc2(1450:1850,1600:1900);
% sitc2 = sitc2(1500:1750,1650:1800);
% figure; imshow(sitc2,[-5 100]);  colormap jet; colorbar;




% 第一天影像sitc1
sitc1(sitc1<15)=0;
sitc1(sitc1>=15)=1;
sitc1(isnan(sitc1))=0;
contour = bwperim(sitc1,8);% 8邻域划分contour
% figure; imshow(contour,[]);
% imwrite(contour,'contour.jpg');
bw1=bwboundaries(sitc1); %获取contour里面有多少个连通体
num=size(bw1,1); %行数，多少个连通体
siw_bw1=zeros(size(sitc1));
indxxx1=0;
for i=1:num %遍历每一个连通体
    data1=bw1{i,1};   %得到每个连通体轮廓线的坐标，一个N*2的矩阵,此坐标为本地图像坐标
    if length(data1)>500 %如果这个连通体的长度大于500，可以视为是整张图像中最主要的边缘线
    indx1=size(sitc1,1).*(data1(:,2)-1)+data1(:,1);%二维数组和一维数组的坐标转换
    siw_bw1(indx1)=1; %标记边缘线的
    indxxx1=[indxxx1,i]; %标记主要的边缘线是哪个连通体中
    end
end
figure;imshow(siw_bw1,[ ],'InitialMagnification','fit');title('data1');

% 第二天影像sitc2
sitc2(sitc2<15)=0;
sitc2(sitc2>=15)=1;
sitc2(isnan(sitc2))=0;
contour2 = bwperim(sitc2,8);% 8邻域划分contour
% figure; imshow(contour2,[]);
bw2=bwboundaries(sitc2); 
num=size(bw2,1); %行数
siw_bw2=nan(size(sitc2));
indxxx2=0;
for i=1:num
    data2=bw2{i,1};   %得到轮廓线的坐标，一个N*2的矩阵,此坐标为本地图像坐标
    if length(data2)>500
        indx2=size(sitc2,1).*(data2(:,2)-1)+data2(:,1);
        siw_bw2(indx2)=1;
        indxxx2=[indxxx2,i];
    end
end
figure;imshow(siw_bw2,[],'InitialMagnification','fit');title('data2');


%===============原始=============================
distmin=cell(length(indxxx2)-1,1);
for i1=2:length(indxxx2)  %sitc2的边缘线,length(indxxx2)=2;→i1=2
    data22=bw2{indxxx2(i1),1};  %sitc2的边缘线坐标，第二张图轮廓线的坐标
    for i2=1:length(data22) %遍历data2的每个点，这个for循环是计算data2(i2)这个点到所有data1点的最小距离
        for j1=2:length(indxxx1)%length(indxxx1)=2→j1=2
            data11=bw1{indxxx1(j1),1};%表示第一张图sitc1轮廓线的坐标
            % 方法1： 计算地理距离
%                 indx=size(siw_temp,1).*(data2(:,2)-1)+data2(:,1);
%                 x=lon2_new(indx);   y=lat2_new(indx);
%                 lon31=lon2_new(data(i2,1),data(i2,2));  lat31=lat2_new(data(i2,1),data(i2,2)); 
%                 tx = txsite('Name','MathWorks','Latitude',lat31,'Longitude',lon31);
%                 rx = rxsite('Name','Fenway Park','Latitude',y,'Longitude',x);
%                 dist1=distance(tx,rx,'euclidean');

            % 方法2： 计算网格的欧式距离
            x1=data11(:,1);  y1=data11(:,2);  %分别获取sitc1轮廓线的x,y坐标
            x2=data22(i2,1); y2=data22(i2,2); %分别获取sitc2轮廓线的x,y坐标
            dist1 = ((x2-x1).^2+(y2-y1).^2).^0.5; %计算data2(i2)这个点到所有data1点的距离
%             lon31=data(i2,1);  lat31=data(i2,2); 
%             dist1=((lon31-x).^2+(lat31-y).^2).^0.5;
            [distmin{i1-1}(i2),I2]=min(dist1);%data2(i2)这个点到所有data1点的距离中最小的那个值，
            %I2能够找到data(i2)一个点到data2(I2)哪个点有着最小距离
            %基于第2条边缘线上的每个像素点计算其与另一条边缘线的最小距离,
        end
    end
end
%==============================================================
%=================根据I2计算I=================================
distmin=cell(length(indxxx2)-1,1);
for i1=2:length(indxxx2)  %sitc2的边缘线,length(indxxx2)=2;→i1=2
    data22=bw2{indxxx2(i1),1};  %sitc2的边缘线坐标，第二张图轮廓线的坐标
    for i2=1:length(data22) %遍历data2的每个点，这个for循环是计算data2(i2)这个点到所有data1点的最小距离
        if i2 == 476
           for j1=2:length(indxxx1)%length(indxxx1)=2→j1=2
                data11=bw1{indxxx1(j1),1};%表示第一张图sitc1轮廓线的坐标
            % 方法1： 计算地理距离
%                 indx=size(siw_temp,1).*(data2(:,2)-1)+data2(:,1);
%                 x=lon2_new(indx);   y=lat2_new(indx);
%                 lon31=lon2_new(data(i2,1),data(i2,2));  lat31=lat2_new(data(i2,1),data(i2,2)); 
%                 tx = txsite('Name','MathWorks','Latitude',lat31,'Longitude',lon31);
%                 rx = rxsite('Name','Fenway Park','Latitude',y,'Longitude',x);
%                 dist1=distance(tx,rx,'euclidean');

            % 方法2： 计算网格的欧式距离
                x1=data11(:,1);  y1=data11(:,2);  %分别获取sitc1轮廓线的x,y坐标
                x2=data22(i2,1); y2=data22(i2,2); %分别获取sitc2轮廓线的x,y坐标
                dist1 = ((x2-x1).^2+(y2-y1).^2).^0.5; %计算data22(i2)这个点到所有data11点的距离
                [distmin{i1-1}(i2),I2]=min(dist1);%data22(i2)这个点到所有data11点的距离中最小的那个值，即data11(I2)
           end
           break
        else
            continue;
        end
    end
end
% ====================================================================


histData=distmin{:}.*3.125; % distmin表示的是网格距离，网格的个数乘以分辨率0.4km
clear('h','histogram');
[maxxx,I]=nanmax(histData); %maxxx表示最大边缘线距离，I表示最大值在矩阵中的序数，可用于追踪最大值在图中的位置
%再利用I回去找I2
minnn=nanmin(histData);%最小边缘线距离
histData(isnan(histData))=[]; %所有基于边缘上的像素点计算的最小距离
nbins=10;
figure;h = histogram(histData,nbins);
hValues=h.Values;
hProsi=hValues./length(histData);%概率直方图
aa=cumsum(hValues);%累积概率直方图
bW=h.BinWidth;
xx=linspace(minnn,maxxx,11);
title('分布直方图');
xlabel('mindistance values(km)'); ylabel('pixel numbers');
xticklabels(xx); 
ax1=gca;  ax1.FontSize=16;   ax1.LineWidth=1.0; ax1.Clipping='on'; ax1.Box='on';
%I代表data2(I)到data1的距离最大
%然后I取值，代入到循环里面，查找I2为多少

%%
figure;imshow(siw_bw2,[ ],'InitialMagnification','fit');  hold on 
rectangle('Position',[88 197 2 2],'FaceColor','b','EdgeColor','b',...   
    'LineWidth',3) %[799 518 10 10] 表示的是点在第799列，518行，矩形大小选10*10个像素  
rectangle('Position',[119 209 2 2],'FaceColor','r','EdgeColor','r',...
    'LineWidth',3)  %[947 712 10 10] 表示的是点在第947列，712行，矩形大小选10*10个像素

figure;imshow(siw_bw1,[ ],'InitialMagnification','fit');  hold on 
rectangle('Position',[88 197 2 2],'FaceColor','b','EdgeColor','b',...
    'LineWidth',3)
rectangle('Position',[119 209 2 2],'FaceColor','r','EdgeColor','r',...
    'LineWidth',3)


% close all %关闭所有图像
%==================================================================================
%%

% 
% % imwrite(fig,'sitc2.jpg');
% 
% 
% % 裁剪图像
% [I,map]=imread('D:/wzx/contour.jpg');
% subplot(121);imshow(I,map);
% %指定剪切区域的大小和位置，剪切，返回xy坐标和裁剪区域
% [x,y,I2,rect]=imcrop(I,map,[500 500 250 250]);%位置和区域大小
% subplot(122);imshow(I2)
% imwrite(I2,'result.jpg');
% 
% 
% latlon_file=[path, 'LongitudeLatitudeGrid-n3125-Arctic3125.hdf'];
% lon2 = double(hdfread(latlon_file,'Longitudes'));
% lat2 = double(hdfread(latlon_file,'Latitudes'));
% lon2_clip=lon2(1000:2500,2000:2400);
% lat2_clip=lat2(1000:2500,2000:2400);
% 
% 
% 
% flag = (lon2>60);
% sitc2(flag)=nan;
% figure;imshow(sitc2,[]);
% 
% %% 指定位置裁剪
% clc
% clear
% [I,map]=imread('D:/wzx/img_edge/af.png');
% figure;
% subplot(121);imshow(I,map);
% %指定剪切区域的大小和位置，剪切，返回xy坐标和裁剪区域
% [x,y,I2,rect]=imcrop(I,map,[1380 1250 1450 1500]);%位置和区域大小
% subplot(122);imshow(I2)
% imwrite(I2,'result.jpg');
% 
% bw1 = bwboundaries(I2);
% num=size(bw1,1); %行数
% siw_bw=zeros(size(I2));
% indxxx1=0;
% for i=1:num
%     data=bw1{i,1};   %得到轮廓线的坐标，一个N*2的矩阵,此坐标为本地图像坐标
%     if length(data)>500
%     indx=size(sitc2_thre,1).*(data(:,2)-1)+data(:,1);
%     siw_bw(indx)=1;
%     indxxx1=[indxxx1,i];
%     end
% end
% figure;imshow(siw_bw,[ ]);
% 
% 
% 
% 
% m_proj('lambert','lon',[-10 20],'lat',[33 48]); 
% 
% 
% minnn=nanmin(sitc2);
% maxxx=nanmax(sitc2);
% 
% lon = lon2((lon2>0)&(lon2<60));
% lat = lat2((lat2>70)&(lon2<80));
% %=========== 读取经纬度hdf文件=============================
% 
% 
% 
% 
% 
% 
% 
% if nanmax(nanmax(lon3))>180 %如果SAR分类结果跨180度，lon2不变
%     lon2=lon2;
% else
%     indx=find(lon2>180);
%     lon2(indx)=lon2(indx)-360;  %如果SAR分类结果不跨180度，lon2改为-180~180
% end
% 
% a =imread('D:\wzx\img_edge\af.png');
% imshow(a)
% bw = im2bw( a );
% contour = bwperim(bw);    
% figure
% imshow(contour)
% 
% ibw = ~bw;
% F1 = imfill(ibw,'holes');
% SE = ones(3); %图像被结构元素SE膨胀
% F2 = imdilate(F1,SE,'same');%膨胀操作
% BW3 = bwperim(F2);
% 
% subplot(1,2,1);imshow(a);title('original iamge');
% subplot(1,2,2);imshow(BW3);title('operated bwperim'); 
% 
% 
% 
% 
% 
% P=hdfread('NISE_SSMISF18_20171124.HDFEOS','Northern Hemisphere','fields','Extent');
% P(P==255)=105;   % Put ocean at top of indices
% P(P>105)=0;
% 
% % According to web site this is is the projection info. I make te radius
% % of my map less than the actual data field though.
% m_proj('azimuthal equal-area','latitude',90,'radius',47,'rectbox','on');
% 
%  clf
% % Plot data as an image
% offs=9036842.762500/6371228; % Convert projection coords to units of earth radii
% image([-offs offs],[offs -offs],P); set(gca,'ydir','normal');
% colormap([.2 .5 .2;  % 0
%            jet(100); % 1-100
%            1 1 1;    % Greenland
%            0 0 0  ;
%            .9 .9 .9; % dry snow
%            .8 .8 .8; % wet snow
%            0 0 .5]); % 105 - now ocean
% caxis([0 105]);
% 
% m_coast('color','k');
% m_grid('linewi',2,'tickdir','out');
% title({'SSM/I Ice cover Nov 24, 2017',''},'fontsize',14,'fontweight','bold');
% 
% hh=colorbar('h');
% set(hh,'tickdir','out');
% xlabel(hh,'% Ice cover');
% 
% 
% figure;
% m_proj('lambert','lon',[-123-25/60 -122-40/60],'lat',[48+42/60 49+9/60]);