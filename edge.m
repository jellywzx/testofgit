clc;
clear all
%=========== ��ȡÿ����ܼ���nc�ļ� =======================
%��ȡ��һ���Ӱ��
path = 'D:\wzx\amsr2\'; 
dstr = '20190315';
ncfile = [path, 'asi-AMSR2-n3125-',dstr,'-v5.4.nc'];
if ~exist(ncfile,'file')
    hdffile = [path, 'asi-AMSR2-n3125-',dstr,'-v5.4.hdf'];
    sitc1 = hdfread(hdffile, '/ASI Ice Concentration');
elseif exist(ncfile,'file')
    sitc = ncread(ncfile, 'z');  
    sitc1=sitc';   % nc�ļ���hdf�ļ���ͼ�Ǻ��������෴�ģ���Ҫת�ã����ܺ;�γ�Ⱦ���ƥ��
end
% figure; 
% fig = imshow(sitc1,[0,100]); colormap jet; colorbar % ��ʾͼ��
sitc1 = sitc1(1450:1850,1600:1900);
% sitc1 = sitc1(1500:1750,1650:1800);
% figure; imshow(sitc1,[-5 100]);  colormap jet; colorbar;

% ========================================================
% ��ȡ�ڶ����Ӱ��
% clc;
% clear all;
path = 'D:\wzx\amsr2\';
dstr2 = '20190323';
hdffile2 = [path, 'asi-AMSR2-n3125-',dstr2,'-v5.4.hdf'];
sitc2 = hdfread(hdffile2, '/ASI Ice Concentration');
% figure; 
% fig = imshow(sitc2,[0,100]); colormap jet; colorbar % ��ʾͼ��
sitc2 = sitc2(1450:1850,1600:1900);
% sitc2 = sitc2(1500:1750,1650:1800);
% figure; imshow(sitc2,[-5 100]);  colormap jet; colorbar;




% ��һ��Ӱ��sitc1
sitc1(sitc1<15)=0;
sitc1(sitc1>=15)=1;
sitc1(isnan(sitc1))=0;
contour = bwperim(sitc1,8);% 8���򻮷�contour
% figure; imshow(contour,[]);
% imwrite(contour,'contour.jpg');
bw1=bwboundaries(sitc1); %��ȡcontour�����ж��ٸ���ͨ��
num=size(bw1,1); %���������ٸ���ͨ��
siw_bw1=zeros(size(sitc1));
indxxx1=0;
for i=1:num %����ÿһ����ͨ��
    data1=bw1{i,1};   %�õ�ÿ����ͨ�������ߵ����꣬һ��N*2�ľ���,������Ϊ����ͼ������
    if length(data1)>500 %��������ͨ��ĳ��ȴ���500��������Ϊ������ͼ��������Ҫ�ı�Ե��
    indx1=size(sitc1,1).*(data1(:,2)-1)+data1(:,1);%��ά�����һά���������ת��
    siw_bw1(indx1)=1; %��Ǳ�Ե�ߵ�
    indxxx1=[indxxx1,i]; %�����Ҫ�ı�Ե�����ĸ���ͨ����
    end
end
figure;imshow(siw_bw1,[ ],'InitialMagnification','fit');title('data1');

% �ڶ���Ӱ��sitc2
sitc2(sitc2<15)=0;
sitc2(sitc2>=15)=1;
sitc2(isnan(sitc2))=0;
contour2 = bwperim(sitc2,8);% 8���򻮷�contour
% figure; imshow(contour2,[]);
bw2=bwboundaries(sitc2); 
num=size(bw2,1); %����
siw_bw2=nan(size(sitc2));
indxxx2=0;
for i=1:num
    data2=bw2{i,1};   %�õ������ߵ����꣬һ��N*2�ľ���,������Ϊ����ͼ������
    if length(data2)>500
        indx2=size(sitc2,1).*(data2(:,2)-1)+data2(:,1);
        siw_bw2(indx2)=1;
        indxxx2=[indxxx2,i];
    end
end
figure;imshow(siw_bw2,[],'InitialMagnification','fit');title('data2');


%===============ԭʼ=============================
distmin=cell(length(indxxx2)-1,1);
for i1=2:length(indxxx2)  %sitc2�ı�Ե��,length(indxxx2)=2;��i1=2
    data22=bw2{indxxx2(i1),1};  %sitc2�ı�Ե�����꣬�ڶ���ͼ�����ߵ�����
    for i2=1:length(data22) %����data2��ÿ���㣬���forѭ���Ǽ���data2(i2)����㵽����data1�����С����
        for j1=2:length(indxxx1)%length(indxxx1)=2��j1=2
            data11=bw1{indxxx1(j1),1};%��ʾ��һ��ͼsitc1�����ߵ�����
            % ����1�� ����������
%                 indx=size(siw_temp,1).*(data2(:,2)-1)+data2(:,1);
%                 x=lon2_new(indx);   y=lat2_new(indx);
%                 lon31=lon2_new(data(i2,1),data(i2,2));  lat31=lat2_new(data(i2,1),data(i2,2)); 
%                 tx = txsite('Name','MathWorks','Latitude',lat31,'Longitude',lon31);
%                 rx = rxsite('Name','Fenway Park','Latitude',y,'Longitude',x);
%                 dist1=distance(tx,rx,'euclidean');

            % ����2�� ���������ŷʽ����
            x1=data11(:,1);  y1=data11(:,2);  %�ֱ��ȡsitc1�����ߵ�x,y����
            x2=data22(i2,1); y2=data22(i2,2); %�ֱ��ȡsitc2�����ߵ�x,y����
            dist1 = ((x2-x1).^2+(y2-y1).^2).^0.5; %����data2(i2)����㵽����data1��ľ���
%             lon31=data(i2,1);  lat31=data(i2,2); 
%             dist1=((lon31-x).^2+(lat31-y).^2).^0.5;
            [distmin{i1-1}(i2),I2]=min(dist1);%data2(i2)����㵽����data1��ľ�������С���Ǹ�ֵ��
            %I2�ܹ��ҵ�data(i2)һ���㵽data2(I2)�ĸ���������С����
            %���ڵ�2����Ե���ϵ�ÿ�����ص����������һ����Ե�ߵ���С����,
        end
    end
end
%==============================================================
%=================����I2����I=================================
distmin=cell(length(indxxx2)-1,1);
for i1=2:length(indxxx2)  %sitc2�ı�Ե��,length(indxxx2)=2;��i1=2
    data22=bw2{indxxx2(i1),1};  %sitc2�ı�Ե�����꣬�ڶ���ͼ�����ߵ�����
    for i2=1:length(data22) %����data2��ÿ���㣬���forѭ���Ǽ���data2(i2)����㵽����data1�����С����
        if i2 == 476
           for j1=2:length(indxxx1)%length(indxxx1)=2��j1=2
                data11=bw1{indxxx1(j1),1};%��ʾ��һ��ͼsitc1�����ߵ�����
            % ����1�� ����������
%                 indx=size(siw_temp,1).*(data2(:,2)-1)+data2(:,1);
%                 x=lon2_new(indx);   y=lat2_new(indx);
%                 lon31=lon2_new(data(i2,1),data(i2,2));  lat31=lat2_new(data(i2,1),data(i2,2)); 
%                 tx = txsite('Name','MathWorks','Latitude',lat31,'Longitude',lon31);
%                 rx = rxsite('Name','Fenway Park','Latitude',y,'Longitude',x);
%                 dist1=distance(tx,rx,'euclidean');

            % ����2�� ���������ŷʽ����
                x1=data11(:,1);  y1=data11(:,2);  %�ֱ��ȡsitc1�����ߵ�x,y����
                x2=data22(i2,1); y2=data22(i2,2); %�ֱ��ȡsitc2�����ߵ�x,y����
                dist1 = ((x2-x1).^2+(y2-y1).^2).^0.5; %����data22(i2)����㵽����data11��ľ���
                [distmin{i1-1}(i2),I2]=min(dist1);%data22(i2)����㵽����data11��ľ�������С���Ǹ�ֵ����data11(I2)
           end
           break
        else
            continue;
        end
    end
end
% ====================================================================


histData=distmin{:}.*3.125; % distmin��ʾ����������룬����ĸ������Էֱ���0.4km
clear('h','histogram');
[maxxx,I]=nanmax(histData); %maxxx��ʾ����Ե�߾��룬I��ʾ���ֵ�ھ����е�������������׷�����ֵ��ͼ�е�λ��
%������I��ȥ��I2
minnn=nanmin(histData);%��С��Ե�߾���
histData(isnan(histData))=[]; %���л��ڱ�Ե�ϵ����ص�������С����
nbins=10;
figure;h = histogram(histData,nbins);
hValues=h.Values;
hProsi=hValues./length(histData);%����ֱ��ͼ
aa=cumsum(hValues);%�ۻ�����ֱ��ͼ
bW=h.BinWidth;
xx=linspace(minnn,maxxx,11);
title('�ֲ�ֱ��ͼ');
xlabel('mindistance values(km)'); ylabel('pixel numbers');
xticklabels(xx); 
ax1=gca;  ax1.FontSize=16;   ax1.LineWidth=1.0; ax1.Clipping='on'; ax1.Box='on';
%I����data2(I)��data1�ľ������
%Ȼ��Iȡֵ�����뵽ѭ�����棬����I2Ϊ����

%%
figure;imshow(siw_bw2,[ ],'InitialMagnification','fit');  hold on 
rectangle('Position',[88 197 2 2],'FaceColor','b','EdgeColor','b',...   
    'LineWidth',3) %[799 518 10 10] ��ʾ���ǵ��ڵ�799�У�518�У����δ�Сѡ10*10������  
rectangle('Position',[119 209 2 2],'FaceColor','r','EdgeColor','r',...
    'LineWidth',3)  %[947 712 10 10] ��ʾ���ǵ��ڵ�947�У�712�У����δ�Сѡ10*10������

figure;imshow(siw_bw1,[ ],'InitialMagnification','fit');  hold on 
rectangle('Position',[88 197 2 2],'FaceColor','b','EdgeColor','b',...
    'LineWidth',3)
rectangle('Position',[119 209 2 2],'FaceColor','r','EdgeColor','r',...
    'LineWidth',3)


% close all %�ر�����ͼ��
%==================================================================================
%%

% 
% % imwrite(fig,'sitc2.jpg');
% 
% 
% % �ü�ͼ��
% [I,map]=imread('D:/wzx/contour.jpg');
% subplot(121);imshow(I,map);
% %ָ����������Ĵ�С��λ�ã����У�����xy����Ͳü�����
% [x,y,I2,rect]=imcrop(I,map,[500 500 250 250]);%λ�ú������С
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
% %% ָ��λ�òü�
% clc
% clear
% [I,map]=imread('D:/wzx/img_edge/af.png');
% figure;
% subplot(121);imshow(I,map);
% %ָ����������Ĵ�С��λ�ã����У�����xy����Ͳü�����
% [x,y,I2,rect]=imcrop(I,map,[1380 1250 1450 1500]);%λ�ú������С
% subplot(122);imshow(I2)
% imwrite(I2,'result.jpg');
% 
% bw1 = bwboundaries(I2);
% num=size(bw1,1); %����
% siw_bw=zeros(size(I2));
% indxxx1=0;
% for i=1:num
%     data=bw1{i,1};   %�õ������ߵ����꣬һ��N*2�ľ���,������Ϊ����ͼ������
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
% %=========== ��ȡ��γ��hdf�ļ�=============================
% 
% 
% 
% 
% 
% 
% 
% if nanmax(nanmax(lon3))>180 %���SAR��������180�ȣ�lon2����
%     lon2=lon2;
% else
%     indx=find(lon2>180);
%     lon2(indx)=lon2(indx)-360;  %���SAR����������180�ȣ�lon2��Ϊ-180~180
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
% SE = ones(3); %ͼ�񱻽ṹԪ��SE����
% F2 = imdilate(F1,SE,'same');%���Ͳ���
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