clear
clc
%$ i=2;
for i=1:2
%$     data_or=xlsread('New Microsoft Excel Worksheet.xlsx',['sheet' num2str(i)]);
    data_or=xlsread('measurement1 - Copy.xlsx',['hole' num2str(i)]);
    data_or=data_or(:,:); %
    data=data_or;
%$    

    [b,a]=butter(2,0.3,'low');
    data(:,2:end)=filtfilt(b,a,data_or(:,2:end));
    [value,inds]=max(data(:,2:end));
%$     figure; scatter(data(inds,1),value);
    inds_1=ones(2,1);
    inds_1(2)=max(inds);
    inds_1(1)=min(inds);
    figure;plot(data(:,1),data(:,2:end));hold on;
    plot(data(:,1),data(:,2:end));hold on;
    y(1)=0;y(2)=1800;
    diff_wave = data(inds_1(2),1)-data(inds_1(1),1);
    diff_ev = 1240/min(value)-1240/max(value);
    for i=1:2
        x_inds=ones(2,1)*data(inds_1(i),1);
        plot(x_inds,y);hold on;
    end
    title(['sheet' num2str(i) '  meV diff:' num2str(diff_wave)])
end
%$energyshift
clc,clear
%$
V_f = 0.0011;%eV0.0012
E_1 = normrnd(0,0.0012,2000,1); $%$$%$519-520
E_2 = normrnd(0,0.0012,2000,1);
$%$ temp_E1 = zeros(100,100);
$%$ temp_E2 = zeros(1,10000);
counter = 1;
delta=zeros(1,4000000);
for k = 1:length(E_1)
%$     temp_E1(:,1)=ones(100,1)*E_1(k);
    for l = 1:length(E_2)
%$         temp_E2(k*l)=E_2(l);
        E1 = E_1(k);
        E2 = E_2(l);
        E_a = ((E1+E2)-sqrt((E1+E2)^2-4*(E1*E2-V_f^2)))/2;
        E_b = ((E1+E2)+sqrt((E1+E2)^2-4*(E1*E2-V_f^2)))/2;
        E_m = min(E_a, E_b);
        
        delta(counter) = (E1-E_m)*1000;
        counter = counter +1;
    end
end
 
%$
 
figure
histfit(delta,9,'kernel')
%$
xlabel('Energy shift (meV)','FontSize',24);
ylabel('Counts','FontSize',24);
set(gca, 'fontsize',26,'LineWidth',1.5)
mean(delta)
std(delta)
pd = fitdist(delta','kernel')
%$
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
data = [experiment_data]
h = histogram(data,6,'FaceColor','#518071','EdgeColor','#006548','Linewidth',1.5);
hold on
testx = 0:0.01:4;
pdf_kernel = pdf(pd,testx);
$%$
plot(testx,pdf_kernel*30,'Color','#C10009','LineWidth',2.5);
xlabel('Energy shift (meV)','FontSize',24);
ylabel('Counts','FontSize',24);
set(gca, 'fontsize',26,'LineWidth',1.5)

%$ Plot trapping data
close all;
clc;
 
fid=fopen('2.txt');
v=textscan(fid,'%.10f','headerlines',15);
fclose(fid);
window = 600;
overlap = window*0.9;
padding = 0;
 
y = v{1};
%$ index=find(y==0);
%$ y(index)=[];
x = 1:length(y);
 
%$ order = 1;
%$ framelen = 11;
%$ sgf = sgolayfilt(y,order,framelen);
 
sd_y = my_sd(y, window, overlap, padding);
sd_x = 1:length(sd_y);
length(x);
length(sd_x);
amp = length(x)/length(sd_x);
%$
figure
plot (x/100000,y,'r');
%$
figure
plot(sd_x/100000,sd_y)
%$
%$ title('Measured Data','FontSize',18,'FontWeight','bold');
%$ grid on;
xlabel('Time (s)','FontSize',14);
ylabel('SD','FontSize',14);
set(gca, 'fontsize',18,'LineWidth',1.5)
%$ hold on;
%$ plot (x/100000,sgf);
 
trap = y(2100000:2200000);
figure
notrap = y(5600000:6600000);
acf = autocorr(trap,'NumLags',1000);
acf_n = autocorr(notrap,'NumLags',1000);
trap_x = ((1:length(acf))./100)'
ntrap_x = ((1:length(acf_n))./100)'
plot((1:length(acf))./100,acf,'LineWidth',3)
xlabel("Lag (ms)",'FontSize',14)
ylabel("Autocorrelation",'FontSize',14)
ylim([0 1])
hold on;
plot(((1:length(acf_n))./100)',acf_n,'LineWidth',3)
hold on;
trap_fit = fit(trap_x,acf,'exp2');
ntrap_fit = fit(ntrap_x,acf_n,'exp2');
%$ plot(trap_x,trap_fit(trap_x),'LineWidth',3);
plot(ntrap_x,ntrap_fit(ntrap_x),'LineWidth',3)
set(gca, 'fontsize',18)
legend('Trapping','No Trapping',"Trapping exp2-fit", "No Trapping exp2-fit")
lgd = legend;
lgd.FontSize = 14;
t4=1000/-(ntrap_fit.d);
t3=1000/-(ntrap_fit.b);
tau = -1/min(abs(t3),abs(t4))

function y = my_sd(signal, windowlength, overlap, zeropad)
delta = windowlength - overlap;
indices = 1:delta:length(signal);
if length(signal) - indices(end) + 1 < windowlength
    if zeropad
        signal(end+1:indices(end)+windowlength-1) = 0;
    else
        indices = indices(1:find(indices+windowlength-1 <= length(signal), 1, 'last'));
    end
end
y = zeros(1, length(indices));
signal = signal.^2;
index = 0;
for i = indices
    index = index+1;
    y(index) = std(signal(i:i+windowlength-1))/(mean(signal(i:i+windowlength-1)));
end

%$ Schrodinger_3D_FEM
%$ This part was edited by Brett Henderson and Hao Zhang
function [bounds] = assign_faces(model, x_ext, y_ext, z_ext)

mesh = model.Mesh;
faces = model.Geometry.NumFaces;  
bounds = zeros(1, 6);
nf = 0;
 
for face = 1:faces
    %$ get the ids of the nodes on a face
    Nf = findNodes(mesh,'region','Face',face);
 
    %$ get the lowest and highest nodes on the face to determine whether the
    %$ face is a bot, top, side, or sphere
    maxz = max(abs(mesh.Nodes(3,Nf)));
    maxx = max(abs(mesh.Nodes(1,Nf)));
    maxy = max(abs(mesh.Nodes(2,Nf)));
    
    eps = 1e-1;
    
    if ((maxz > (z_ext + eps)) || (maxy > (y_ext + eps)) || (maxx > (x_ext + eps)))  % top
        bounds(nf+1) = face;
        nf = nf + 1;
    end
end
clear all;
close all;
clc;
model = createpde();
 
%$ Define the Computational Domain and visualize it
m_star = 0.15;
V_in = 0;    %$ in eV
V_out = 4;  %$ in eV
V_lig = 2;
 
box_width = 12;
bx = box_width;  %$ in nm
by = box_width;  %$ in nm
bz = box_width;  %$ in nm
lig_width = 1;  %$ in nm
 
max_mesh_size = bx / 6;  %$ nm
eval_int = [0, 0.4];
 
cellx_mult = 2.25;  %$ total simulation cell $%$dimension relative to box
celly_mult = 2.25;
cellz_mult = 2.25;
 
cellx = cellx_mult * bx;
celly = celly_mult * by;
cellz = cellz_mult * bz;
 
xx = [linspace(-floor(cellx/2), -(bx/2+lig_width), 5), ...
      linspace(-(bx/2+lig_width), -bx/2, 5), ...
      linspace(-bx/2, bx/2, 5), ...
      linspace(bx/2, bx/2+lig_width, 5), ...
      linspace(bx/2+lig_width, floor(cellx / 2), 5), ...
     ];
 
yy = [linspace(-floor(celly/2), -(by/2+lig_width), 5), ...
      linspace(-(by/2+lig_width), -by/2, 5), ...
      linspace(-by/2, by/2, 5), ...
      linspace(by/2, by/2+lig_width, 5), ...
      linspace(by/2+lig_width, floor(celly / 2), 5), ...
     ];
 
zz = [linspace(-floor(cellz/2), -(bz/2+lig_width), 5), ...
      linspace(-(bz/2+lig_width), -bz/2, 5), ...
      linspace(-bz/2, bz/2, 5), ...
      linspace(bz/2, bz/2+lig_width, 5), ...
      linspace(bz/2+lig_width, floor(cellz / 2), 5), ...
     ];
 
xx = unique(xx);
yy = unique(yy);
zz = unique(zz);
 
[xg, yg, zg] = meshgrid(xx, yy, zz);
                    
Pcube = [xg(:) yg(:), zg(:)];
                
shpCube = alphaShape(Pcube); 
                           
 figure
 plot(shpCube,'FaceAlpha',0.4)
 title('alphaShape: Cube')
 
%$ Triangulate the Domain
 
[tri,loc] = alphaTriangulation(shpCube);
[gCube,msh] = geometryFromMesh(model,loc',tri');
 
figure
pdegplot(model,'FaceAlpha',0.5,'CellLabels','on')
title('PDEModel: Cube')
 
%$ Assign the Inner Cube to a Second Domain
nodes = msh.Nodes;
elements = msh.Elements;
 
elemXCoords = reshape(nodes(1,elements),4,[]);
elemXCoordsGeometricCenter = mean(elemXCoords);
elemYCoords = reshape(nodes(2,elements),4,[]);
elemYCoordsGeometricCenter = mean(elemYCoords);
elemZCoords = reshape(nodes(3,elements),4,[]);
elemZCoordsGeometricCenter = mean(elemZCoords);
 
ElementIdToRegionId = ones(1,size(elements,2));
 
% ligand + cube
idx = (abs(elemXCoordsGeometricCenter) <= bx / 2 + lig_width) ...
      & (abs(elemYCoordsGeometricCenter) <= by / 2 + lig_width) ...
      & (abs(elemZCoordsGeometricCenter) <= bz / 2 + lig_width);
ElementIdToRegionId(idx) = 3;
 
% cube
idx = (abs(elemXCoordsGeometricCenter) <= bx / 2) ...
      & (abs(elemYCoordsGeometricCenter) <= by / 2) ...
      & (abs(elemZCoordsGeometricCenter) <= bz / 2);
ElementIdToRegionId(idx) = 2;
 
nm_to_bohr = 18.8973; 
lig_width_bohr = lig_width * nm_to_bohr;
bx_bohr = bx * nm_to_bohr;
by_bohr = by * nm_to_bohr;
bz_bohr = bz * nm_to_bohr;
 
model2 = createpde;
geometryFromMesh(model2,nodes*nm_to_bohr,elements,ElementIdToRegionId);
pdegplot(model2,'CellLabels','on', 'FaceLabels','on', 'FaceAlpha',0.5);
 
bounds = assign_faces(model2, bx_bohr/2+lig_width_bohr, by_bohr/2+lig_width_bohr, bz_bohr/2+lig_width_bohr);
applyBoundaryCondition(model2,'dirichlet','face',bounds,'u',0);
 

h_bar = 1;
Eh_to_eV = 27.211386245988;
c = 1/(2 * m_star);
a_in = V_in / Eh_to_eV;
a_lig = V_lig / Eh_to_eV;
a_out = V_out / Eh_to_eV;
 
specifyCoefficients(model2,'m',0,'f',0,'a',a_out,'d',1,'c',c,'Cell',1);  $%$ boundaries
specifyCoefficients(model2,'m',0,'f',0,'a',a_in,'d',1,'c',c,'Cell',2);  $%$ well
specifyCoefficients(model2,'m',0,'f',0,'a',a_lig,'d',1,'c',c,'Cell',3);   $%$ ligand

 
generateMesh(model2,'Hmax',max_mesh_size * nm_to_bohr);
elementIDsCell1 = findElements(model2.Mesh,'region','Cell',1);
elementIDsCell2 = findElements(model2.Mesh,'region','Cell',2);
 
figure
pdemesh(model2.Mesh.Nodes, ...
        model2.Mesh.Elements(:,elementIDsCell1), ...
        'FaceColor','red', 'FaceAlpha', 0.5)
hold on
pdemesh(model2.Mesh.Nodes, ...
        model2.Mesh.Elements(:,elementIDsCell2), ...
        'FaceColor','green')

eval_int = eval_int / Eh_to_eV;
result = solvepdeeig(model2, eval_int);
result.Eigenvalues * Eh_to_eV

xp = -200:10:200; yp = -200:10:200;
[XX, YY] = meshgrid(xp,yp);
ZZ = zeros(size(XX));
 
uintrp = interpolateSolution(result,XX,YY,ZZ,1);
ures = reshape(uintrp,size(XX));
surf(XX,YY,ures)
shading flat
shading interp
xlabel('x'); ylabel('y'); zlabel('z');
title('Ground State')
V = result.Eigenvectors;
subplot(2,2,1)
pdeplot3D(model,'ColorMapData',V(:,1))
title('Ground State')
subplot(2,2,2)
pdeplot3D(model,'ColorMapData',V(:,2))
title('First Excited State')
subplot(2,2,3)
pdeplot3D(model,'ColorMapData',V(:,3))
title('Second Excited State')
