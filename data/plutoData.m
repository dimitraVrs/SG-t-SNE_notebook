% Points for section CLUSTERS' SIZE

muSize=zeros(50,50);
muSize(26:50,1)=100;
sigmaSize=ones(50,50);
sigmaSize(1:25,:)=10;

clSize=normrnd(muSize,sigmaSize,50,50);

figure(2);
scatter(clSize(1:25,1),clSize(1:25,2),[],[0 0.4470 0.7410])
hold on
scatter(clSize(26:50,1),clSize(26:50,2),[],[0.8500 0.3250 0.0980]);

pSize=normpdf(clSize,muSize,sigmaSize);

for i=1:50
    pSize(i,i)=0;
end

nnzSize=nnz(pSize);

colToAddSize=zeros(50,nnzSize-50);
tempSize2=horzcat(pSize,colToAddSize);
rowsCluster1Size=zeros(fix((nnzSize-50)/2),nnzSize);
cl1Size=vertcat(tempSize2(1:25,:),rowsCluster1Size);
rowsCluster2Size=zeros(fix((nnzSize-50)/2)+mod(nnzSize-50,2),nnzSize);
cl2Size=vertcat(tempSize2(26:50,:),rowsCluster2Size);

pointsSizeSp=vertcat(cl1Size,cl2Size);

pointsSizeIn=sparse(pointsSizeSp);

mmwrite("sizeData.mtx",pointsSizeIn);

% Points for section CLUSTERS' DISTANCE

muD=zeros(25,75);
sigmaD=ones(25,75);

s = rng; % save state of rng
cluster1D=normrnd(muD,sigmaD,25,75);
rng(s); % restore state of rng
cluster2D=normrnd(muD,sigmaD,25,75); % cluster2 has the same random numbers as cluster1
rng(s);
cluster3D=normrnd(muD,sigmaD,25,75);

cluster2D(:,1)=cluster2D(:,1)+10;
mu2D(:,1)=muD(:,1)+10;

cluster3D(:,1)=cluster3D(:,1)+40; % cluster 1 & 2 have distance 10
                                % cluster 1 & 3 have distance 40

mu3D(:,1)=muD(:,1)+40;

figure(3);
scatter(cluster1D(:,1),cluster1D(:,2),[],[0 0.4470 0.7410])
hold on
scatter(cluster2D(:,1),cluster2D(:,2),[],[0.8500 0.3250 0.0980]);
hold on
scatter(cluster3D(:,1),cluster3D(:,2),[],[0.4660 0.6740 0.1880]);

pdf1D=normpdf(cluster1D,muD,sigmaD);
pdf2D=normpdf(cluster2D,mu2D,sigmaD);
pdf3D=normpdf(cluster3D,mu3D,sigmaD);

tempD1=vertcat(pdf1D,pdf2D);
pointsD=vertcat(tempD1,pdf3D);

for i=1:75
    pointsD(i,i)=0;
end

nnzD=nnz(pointsD);

mmwrite("distDataTest.mtx",pointsD);

colToAddD=zeros(75,nnzD-75);
tempD2=horzcat(pointsD,colToAddD);
rowsCluster1D=zeros(fix((nnzD-75)/3),nnzD);
cl1D=vertcat(tempD2(1:25,:),rowsCluster1D);
rowsCluster2D=zeros(fix((nnzD-75)/3),nnzD);
cl2D=vertcat(tempD2(26:50,:),rowsCluster2D);
rowsCluster3D=zeros(fix((nnzD-75)/3)+mod(nnzD-75,3),nnzD);
cl3D=vertcat(tempD2(51:75,:),rowsCluster3D);

tempD3=vertcat(cl1D,cl2D);
pointsD=vertcat(tempD3,cl3D);

pointsDIn=sparse(pointsD);

mmwrite("distData.mtx",pointsDIn);

% Points for section TOPOLOGY

% Circles

muTop=zeros(25,50);
sigmaTop=ones(25,50);

s = rng;
cluster1Top=normrnd(muTop,sigmaTop,25,50);
rng(s);
cluster2Top=normrnd(muTop,sigmaTop*50,25,50);

figure(1);
scatter(cluster1Top(:,1),cluster1Top(:,2),[],[0 0.4470 0.7410])
hold on
scatter(cluster2Top(:,1),cluster2Top(:,2),[],[0.8500 0.3250 0.0980]);

pdf1Top=normpdf(cluster1Top,muTop,sigmaTop);
pdf2Top=normpdf(cluster2Top,muTop,sigmaTop*50);

pointsTop=vertcat(pdf1Top,pdf2Top);

for i=1:50
    pointsTop(i,i)=0;
end

nnzTop=nnz(pointsTop);

colToAddTop=zeros(50,nnzTop-50);
tempTop2=horzcat(pointsTop,colToAddTop);
rowsCluster1Top=zeros(fix((nnzTop-50)/2),nnzTop);
cl1Top=vertcat(tempTop2(1:25,:),rowsCluster1Top);
rowsCluster2Top=zeros(fix((nnzTop-50)/2)+mod(nnzTop-50,2),nnzTop);
cl2Top=vertcat(tempTop2(26:50,:),rowsCluster2Top);

% tempcl1Top=6*ones(size(cl1Top,1)-25,50);
% tempcl1Top(1:fix((size(cl1Top,1)-25)/2),:)=-tempcl1Top(1:fix((size(cl1Top,1)-25)/2),:);
%
% tempcl2Top=6*50*ones(size(cl2Top,1)-25,50);
% tempcl2Top(1:fix((size(cl2Top,1)-25)/2)+mod(size(cl2Top,1)-25,2),:)=-tempcl2Top(1:fix((size(cl2Top,1)-25)/2),:);
%
% temp2cl1Top=vertcat(cluster1Top,tempcl1Top);
% xcl1Top=horzcat(temp2cl1Top,zeros(size(temp2cl1Top,1),nnzTop-50));
%
% temp2cl2Top=vertcat(cluster2Top,tempcl2Top);
% xcl2Top=horzcat(temp2cl2Top,zeros(size(temp2cl2Top,1),nnzTop-50));

pointsTopSp=vertcat(cl1Top,cl2Top);

pointsTopIn=sparse(pointsTopSp);

mmwrite("topData.mtx",pointsTopIn);

% Trefoil knot

% muTop=zeros(25,50);
% sigmaTop=ones(25,50);
% 
% s = rng;
% cluster1Top=normrnd(muTop,sigmaTop,25,50);
% rng(s);
% cluster2Top=normrnd(muTop,sigmaTop*50,25,50);

clusterTr=zeros(50,50);
cTr=zeros(50,1);
for i=1:50
    t=i*2*pi/50;
    clusterTr(i,1)=sin(t)+2*sin(2*t);
    clusterTr(i,2)=cos(t)-2*cos(2*t);
    clusterTr(i,3)=-sin(3*t);
    cTr(i,1)=t;
end

figure(7);
scatter3(clusterTr(:,1),clusterTr(:,2),clusterTr(:,3),[],cTr)

muTr=zeros(50,50);
sigmaTr=ones(50,50);
pdfTr=normpdf(clusterTr,muTr,sigmaTr);

for i=1:50
    pdfTr(i,i)=0;
end

nnzTr=nnz(pdfTr);

colToAddTr=zeros(50,nnzTr-50);
tempTr2=horzcat(pdfTr,colToAddTr);
rowsTr=zeros(nnzTr-50,nnzTr);
pointsTr=vertcat(tempTr2,rowsTr);

pointsTrIn=sparse(pointsTr);

mmwrite("topData2.mtx",pointsTrIn);

% Points for section NOISE

muN=zeros(75,75);
sigmaN=ones(75,75);

clusterN=normrnd(muN,sigmaN,75,75);

figure(4);
scatter(clusterN(:,1),clusterN(:,2),[],[0 0.4470 0.7410])

pN=normpdf(clusterN,muN,sigmaN);

for i=1:75
    pN(i,i)=0;
end

nnzN=nnz(pN);

colToAddN=zeros(75,nnzN-75);
tempN2=horzcat(pN,colToAddN);
rowsToAddN=zeros(nnzN-75,nnzN);
pointsN=vertcat(tempN2,rowsToAddN);

pointsNIn=sparse(pointsN);

mmwrite("noiseData.mtx",pointsNIn);

% Points for section SHAPES

% Ellipse

muSh=zeros(50,50);
sigmaSh=ones(50,50);


clusterSh=normrnd(muSh,sigmaSh,50,50);

for i=1:50
    for j=1:50
        clusterSh(i,j)=clusterSh(i,j)/j;
        muSh(i,j)=muSh(i,j)/j;
        sigmaSh(i,j)=sigmaSh(i,j)/j;
    end
end

figure(5);
scatter(clusterSh(:,1),clusterSh(:,2),[],[0 0.4470 0.7410])

pSh=normpdf(clusterSh,muSh,sigmaSh);

for i=1:50
    pSh(i,i)=0;
end

nnzSh=nnz(pSh);

colToAddSh=zeros(50,nnzSh-50);
tempSh2=horzcat(pSh,colToAddSh);
rowsToAddSh=zeros(nnzSh-50,nnzSh);
pointsSh=vertcat(tempSh2,rowsToAddSh);

pointsShIn=sparse(pointsSh);

mmwrite("shapesData1.mtx",pointsShIn);

% Parallel lines

m2Sh=zeros(25,50);
sigma2Sh=ones(25,50);


cluster1Sh=normrnd(m2Sh,sigma2Sh,25,50);
cluster2Sh=normrnd(m2Sh,sigma2Sh,25,50);

s=0.03*50;

for i=1:25
    cluster1Sh(i,1)=s*cluster1Sh(i,1)+i;
    cluster2Sh(i,1)=s*cluster2Sh(i,1)+i+50/5;
    for j=2:50
        cluster1Sh(i,j)=s*cluster1Sh(i,j)+i;
        cluster2Sh(i,j)=s*cluster2Sh(i,j)+i-50/5;
    end
end

figure(6);
scatter(cluster1Sh(:,1),cluster1Sh(:,2),[],[0 0.4470 0.7410])
hold on
scatter(cluster2Sh(:,1),cluster2Sh(:,2),[],[0.8500 0.3250 0.0980]);

mu1Sh=zeros(25,50);
mu2Sh=zeros(25,50);
for i=1:25
    mu1Sh(i,1)=s*muSh(i,1)+i;
    mu2Sh(i,1)=s*muSh(i,1)+i+50/5;
    for j=2:50
        mu1Sh(i,j)=s*muSh(i,j)+i;
        mu2Sh(i,j)=s*muSh(i,j)+i-50/5;
    end
end
p1Sh=normpdf(cluster1Sh,mu1Sh,sigma2Sh*s);
p2Sh=normpdf(cluster2Sh,mu2Sh,sigma2Sh*s);

points2Sh=vertcat(p1Sh,p2Sh);

for i=1:50
    points2Sh(i,i)=0;
end

nnz2Sh=nnz(points2Sh);

colToAdd2Sh=zeros(50,nnz2Sh-50);
tempSh2b=horzcat(points2Sh,colToAdd2Sh);
rowsCluster1Sh=zeros(fix((nnz2Sh-50)/2),nnz2Sh);
cl1Sh=vertcat(tempSh2b(1:25,:),rowsCluster1Sh);
rowsCluster2Sh=zeros(fix((nnz2Sh-50)/2)+mod(nnz2Sh-50,2),nnz2Sh);
cl2Sh=vertcat(tempSh2b(26:50,:),rowsCluster2Sh);

points2ShSp=vertcat(cl1Sh,cl2Sh);

points2ShIn=sparse(points2ShSp);

mmwrite("shapesData2.mtx",points2ShIn);

% 2 clusters with the same size

muSize=zeros(50,50);
muSize(26:50,1)=10;
sigmaSize=ones(50,50);

clSize=normrnd(muSize,sigmaSize,50,50);

figure(11);
scatter(clSize(1:25,1),clSize(1:25,2),[],[0 0.4470 0.7410])
hold on
scatter(clSize(26:50,1),clSize(26:50,2),[],[0.8500 0.3250 0.0980]);
title('2 clusters with the same size')

pSize=normpdf(clSize,muSize,sigmaSize);

for i=1:50
    pSize(i,i)=0;
end

nnzSize=nnz(pSize);

pointsSizeSp= makeSparse2Cl(pSize,nnzSize)

pointsSizeIn=sparse(pointsSizeSp);

mmwrite("same_size.mtx",pointsSizeIn);

% points evenly in circle

clusterCircle=zeros(50,50);
c=zeros(50,1);
for i=1:50
        t=2*pi*i/50;
        clusterCircle(i,1)=cos(t);
        clusterCircle(i,2)=sin(t);
        c(i,1)=t;
end

figure(12);
scatter(clusterCircle(:,1),clusterCircle(:,2),[],c)
title('points evenly in circle')

muCircle=mean(clusterCircle,1);

sigmaCircle=std(clusterCircle,1);

pCircle=normpdf(clusterCircle,muCircle,sigmaCircle);

for i=1:50
    pCircle(i,i)=0;
end

nnzCircle=nnz(pCircle);

pCircleSp= makeSparse1Cl(pCircle,nnzCircle);

pointsCircleIn=sparse(pCircleSp);

mmwrite("circle.mtx",pointsCircleIn);

% points randomly in circle

clusterCircle=zeros(50,50);
c=zeros(50,1);
for i=1:50
        t=2*pi*rand();
        clusterCircle(i,1)=cos(t);
        clusterCircle(i,2)=sin(t);
        c(i,1)=t;
end

figure(13);
scatter(clusterCircle(:,1),clusterCircle(:,2),[],c)
title('points randomly in circle')

muCircle=mean(clusterCircle,1);
sigmaCircle=std(clusterCircle,1);

pCircle2=normpdf(clusterCircle,muCircle,sigmaCircle);

for i=1:50
    pCircle2(i,i)=0;
end

nnzCircle2=nnz(pCircle2);

pCircleSp2= makeSparse1Cl(pCircle2,nnzCircle2);

pointsCircleIn=sparse(pCircleSp2);

mmwrite("random-circle.mtx",pointsCircleIn);