%% Non-convex mesh, eipsilon = 1
n = 4; eipsilon = 1; error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1); H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
data1 = load('P1P0P0_FracturedDarcy_Polygon_100.mat'); node1 = data1.Node; elem1 = data1.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(1), errorb(1), H1p(1), H1pGamma(1), maxpGamma(1), dof(1)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data2 = load('P1P0P0_FracturedDarcy_Polygon_400.mat'); node1 = data2.Node; elem1 = data2.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(2), errorb(2), H1p(2), H1pGamma(2), maxpGamma(2), dof(2)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data3 = load('P1P0P0_FracturedDarcy_Polygon_1600.mat'); node1 = data3.Node; elem1 = data3.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(3), errorb(3), H1p(3), H1pGamma(3), maxpGamma(3), dof(3)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data4 = load('P1P0P0_FracturedDarcy_Polygon_6400.mat'); node1 = data4.Node; elem1 = data4.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(4), errorb(4), H1p(4), H1pGamma(4), maxpGamma(4), dof(4)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1); H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 1: CVT mesh when eipsilon is 1');
colname = {'#DOF',         'error0', 'rate'  'errorb', 'rate', 'H1p', 'rate',  'H1pGamma', 'rate', 'maxpGamma','rate'};
disptable(colname,dof,[], error0,'%0.4e', error0Rate, [], errorb,'%0.4e', errorbRate, [], ...
            H1p, '%0.4e', H1pRate, [], H1pGamma,'%0.4e', H1pGammaRate, [], maxpGamma,'%0.4e', maxpGammaRate, []);
clear all

%% Non-convex mesh, eipsilon = 1000
n = 4; eipsilon = 1000; error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1); H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
data1 = load('P1P0P0_FracturedDarcy_Polygon_100.mat'); node1 = data1.Node; elem1 = data1.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(1), errorb(1), H1p(1), H1pGamma(1), maxpGamma(1), dof(1)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data2 = load('P1P0P0_FracturedDarcy_Polygon_400.mat'); node1 = data2.Node; elem1 = data2.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(2), errorb(2), H1p(2), H1pGamma(2), maxpGamma(2), dof(2)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data3 = load('P1P0P0_FracturedDarcy_Polygon_1600.mat'); node1 = data3.Node; elem1 = data3.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(3), errorb(3), H1p(3), H1pGamma(3), maxpGamma(3), dof(3)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data4 = load('P1P0P0_FracturedDarcy_Polygon_6400.mat'); node1 = data4.Node; elem1 = data4.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(4), errorb(4), H1p(4), H1pGamma(4), maxpGamma(4), dof(4)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1); H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 2: CVT mesh when eipsilon is 1000');
colname = {'#DOF',         'error0', 'rate'  'errorb', 'rate', 'H1p', 'rate',  'H1pGamma', 'rate', 'maxpGamma','rate'};
disptable(colname,dof,[], error0,'%0.4e', error0Rate, [], errorb,'%0.4e', errorbRate, [], ...
            H1p, '%0.4e', H1pRate, [], H1pGamma,'%0.4e', H1pGammaRate, [], maxpGamma,'%0.4e', maxpGammaRate, []);
clear all

%% Non-convex mesh, eipsilon = 0.001
n = 4; eipsilon = 0.001; error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1); H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
data1 = load('P1P0P0_FracturedDarcy_Polygon_100.mat'); node1 = data1.Node; elem1 = data1.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(1), errorb(1), H1p(1), H1pGamma(1), maxpGamma(1), dof(1)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data2 = load('P1P0P0_FracturedDarcy_Polygon_400.mat'); node1 = data2.Node; elem1 = data2.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(2), errorb(2), H1p(2), H1pGamma(2), maxpGamma(2), dof(2)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data3 = load('P1P0P0_FracturedDarcy_Polygon_1600.mat'); node1 = data3.Node; elem1 = data3.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(3), errorb(3), H1p(3), H1pGamma(3), maxpGamma(3), dof(3)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

data4 = load('P1P0P0_FracturedDarcy_Polygon_6400.mat'); node1 = data4.Node; elem1 = data4.Elem;
node2 = [node1(:,1) -node1(:,2)]; elem2 = cell(length(elem1),1);
for i = 1:length(elem1), elem2{i} = [elem1{i}(1) fliplr(elem1{i}(2:end))]; end
[error0(4), errorb(4), H1p(4), H1pGamma(4), maxpGamma(4), dof(4)] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon);

error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1); H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 3: CVT mesh when eipsilon is 0.001');
colname = {'#DOF',         'error0', 'rate'  'errorb', 'rate', 'H1p', 'rate',  'H1pGamma', 'rate', 'maxpGamma','rate'};
disptable(colname,dof,[], error0,'%0.4e', error0Rate, [], errorb,'%0.4e', errorbRate, [], ...
            H1p, '%0.4e', H1pRate, [], H1pGamma,'%0.4e', H1pGammaRate, [], maxpGamma,'%0.4e', maxpGammaRate, []);
clear all