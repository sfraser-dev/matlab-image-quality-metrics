%%% Implements the quality evaluation index for image restoration
%%% described by Jimenez et al (1-d version). 
%%% All 4 structures, noiseVar=0.0269;
%%% Blocks only, noiseVar=0.021
%%% Triangles only, noiseVar=0.039
%%% Curves only, noiseVar=0.029
%%% Needles only, noiseVar=0.06

function simArrNoisy=jimenez14spBLOCKS(noiseVar,filterNameStr,ws_thr,par2,par3);

base=10; top=30; skip=8;
simArr=[];
simArr=[simArr ones(1,skip)*base];
simArr=[simArr ones(1,42)*top];
simArr=[simArr ones(1,skip)*base];
triGrad=0.6;
simArr=[simArr base:triGrad:top];
simArr=[simArr top-1:-1*triGrad:base];
simArr=[simArr ones(1,skip)*base];
x_array=1:26; dampVal=-0.2;
simArr=[simArr base+(top-base)*(1-exp(dampVal*x_array))];
x_array=reverse(x_array);
simArr=[simArr base+(top-base)*(1-exp(dampVal*x_array))];
clear x_array;
simArr=[simArr ones(1,skip)*base];
needleSep=2; topLength=2; theGradient=5;
Ntop=25;	% needle top.
needle=[ones(1,needleSep)*base ... 
	base:theGradient:Ntop ...
	ones(1,topLength)*Ntop ...
  	Ntop:-1*theGradient:base ...
	ones(1,needleSep)*base]; 
simArr=[simArr needle needle needle needle];
simArr=[simArr ones(1,skip)*base];

%simArr=simArr(1:56); % 1 block
%simArr=simArr(57:128); % 1 triangle
%simArr=simArr(129:188); % 1 curve
%simArr=simArr(189:256);	% 1 set of needles
% 4 blocks
simArr=[simArr(1:56) ones(1,8)*base simArr(1:56) ones(1,8)*base simArr(1:56) ones(1,8)*base simArr(1:56) ones(1,8)*base]; 
% 3 triangles
%simArr=[ones(1,4)*base simArr(57:128) ones(1,12)*base simArr(57:128) ones(1,12)*base simArr(57:128) ones(1,12)*base];
% 3 curves
%simArr=[ones(1,19)*base simArr(129:188) ones(1,19)*base simArr(129:188) ones(1,19)*base simArr(129:188) ones(1,19)*base];
% 3 sets of needles
%simArr=[ones(1,9)*base simArr(189:256) ones(1,14)*base simArr(189:256) ones(1,14)*base simArr(189:256) ones(1,14)*base];
fprintf('length simArr = %f\n',length(simArr));

%%% Add noise to simulated blocky image.
%randn('seed',8888); rand('seed',8888);
simArrNoisy=simArr + sqrt(noiseVar)*simArr.*(rand(size(simArr))-0.5); % speckle
% Calculate the standard deviation of the 1-d profile.
the_std=sqrt(sum(((simArrNoisy-simArr).^2))/prod(size(simArrNoisy)));

%%%
%ws_thr=13; filterType='median';			% median 13
%ws_thr=17; par2=7; par3=15; filterType='elee';		% elee 9 20 10 (g0.62)
%ws_thr=0.11; filterType='daub4';			% daub4 ht 0.11
%ws_thr=95; par2=0.725;	filterType='yu at';		% yu at ht th0+ 95 0.725
%ws_thr=95; par2=1.3;filterType='yu ds';		% yu ds st th0+ 95 1.3
%ws_thr=0.13; filterType='cs';				% cs daub4 ht 0.13 

%%% Remove noise filtering.
	%filtered=AV_M(simArrNoisy,ws_thr);
filtered=MED_M(simArrNoisy,ws_thr);
%simArrNoisy=simArrNoisy+1;filtered=ELEE_M_fl(simArrNoisy,ws_thr,par2,par3);simArrNoisy=simArrNoisy-1;filtered=filtered-1;
	%filtered=wavfilt_ST_PO_log(simArrNoisy,ws_thr);
%filtered=wavfilt_HT_PO_log(simArrNoisy,ws_thr);
    	%filtered=wavfilt_atrousSTlog(simArrNoisy,ws_thr);
    %filtered=wavfilt_atrousHTlog(simArrNoisy,ws_thr);
	%filtered=wavfilt_atrousSSHHTlog(simArrNoisy,ws_thr);
	%filtered=wavfilt_atrousSHHHTlog(simArrNoisy,ws_thr);
	%filtered=datasieve_STlog(simArrNoisy,ws_thr);
	%filtered=datasieve_HTlog(simArrNoisy,ws_thr);
	%filtered=datasieve_SSHTlog(simArrNoisy,ws_thr);
	%filtered=datasieve_SHHHTlog(simArrNoisy,ws_thr);
%filtered=yu_at4SPEK_1d(simArrNoisy,0.001,ws_thr,0.05,par2,'HT','th0+');
%filtered=yu_datasieve4SPEK_1d(simArrNoisy,0.01,ws_thr,1,par2,'ST','th0+');
%filtered=cycleSpin(simArrNoisy,ws_thr);

%%% Print std of noise to the screen.
fprintf('STD of profile: %f\n',the_std);

%%% OLD STYLE ES + NR for blocks only.
ht_sa=simArr(1,20:40);
ht_san=simArrNoisy(1,20:40);
ht_f=filtered(1,20:40);
%plot(ht_sa); hold on; plot(ht_san,'r'); plot(ht_f); hold off;
fprintf('std ht_sa  = %f\n',std(ht_sa));
fprintf('std ht_san = %f\n',std(ht_san));
fprintf('std ht_f   = %f\n',std(ht_f)); 
fprintf('speck reduction percentage: %f\n',100*(1-(std(ht_f)/std(ht_san))));
et_sa=simArr(1,48:53)
et_san=simArrNoisy(1,48:53)
et_f=filtered(1,48:53)
es_sa=mean(et_sa(1,1:3)) / mean(et_sa(4:6));
es_san=mean(et_san(1,1:3)) / mean(et_san(4:6));
es_f=mean(et_f(1,1:3)) / mean(et_f(4:6));
fprintf('es_sa  = %f\n',es_sa);
fprintf('es_san = %f\n',es_san);
fprintf('es_f   = %f\n',es_f);
fprintf('edge integrity percentage: %f\n',100*(es_f/es_san));
ne_sa= sqrt( var(et_sa(1,1:3)) + var(et_sa(1,4:6)) );
ne_san= sqrt( var(et_san(1,1:3)) + var(et_san(1,4:6)) );
ne_f= sqrt( var(et_f(1,1:3)) + var(et_f(1,4:6)) );
fprintf('ne_sa  = %f\n',ne_sa);
fprintf('ne_san = %f\n',ne_san);
fprintf('ne_f   = %f\n',ne_f);
fprintf('speckle red abt edges: %f\n',100*(1-ne_f/ne_san));

%%% Calculate the Smoothing Index (SI).
yRugStore=dAlphaY(simArr,filtered);
rugosity=sum(yRugStore)/prod(size(simArr));
SI=exp(-1*rugosity);
fprintf('SI:  %f\n',SI);

%%% Calculate the Fidelity Index (FI).
FI=fidelity(simArr,filtered);
fprintf('FI:  %f\n',FI);

%%% Global Restoration Index (GRI)
GRI=sqrt(FI*SI);
fprintf('GRI: %f\n',GRI);

%%% PSNR
psnrNoisy=psnr(simArrNoisy,simArr);
psnrFiltered=psnr(filtered,simArr);
fprintf('noisy psnr = %f\n',psnrNoisy);
fprintf('filtered psnr = %f\n',psnrFiltered);
fprintf('psnr improvement = %f\n',psnrFiltered-psnrNoisy);

%%% Image/Signal Fidelity Metric
isfmNoisy=imgSigFidMet(simArrNoisy,simArr);
isfmFiltered=imgSigFidMet(filtered,simArr);
fprintf('noisy img/sig fid met = %f\n',isfmNoisy);
fprintf('filtered img/sig fid met = %f\n',isfmFiltered);
fprintf('img/sig fid met improvement = %f\n\n',isfmFiltered-isfmNoisy);

%%% Outputting data to file.
fid=fopen('jimenez14spBLOCKS.dat','a');
fprintf(fid,'%s\n',filterNameStr);
fprintf(fid,'noiseVar = %f\n',noiseVar);
fprintf(fid,'Noise STD: %f\n',the_std);
fprintf(fid,'std ht_sa  = %f\n',std(ht_sa));
fprintf(fid,'std ht_san = %f\n',std(ht_san));
fprintf(fid,'std ht_f   = %f\n',std(ht_f)); 
fprintf(fid,'speck reduction percentage: %f\n',100*(1-(std(ht_f)/std(ht_san))));
fprintf(fid,'es_sa  = %f\n',es_sa);
fprintf(fid,'es_san = %f\n',es_san);
fprintf(fid,'es_f   = %f\n',es_f);
fprintf(fid,'edge integrity percentage: %f\n',100*(es_f/es_san));
fprintf(fid,'ne_sa  = %f\n',ne_sa);
fprintf(fid,'ne_san = %f\n',ne_san);
fprintf(fid,'ne_f   = %f\n',ne_f);
fprintf(fid,'speckle red abt edges: %f\n',100*(1-ne_f/ne_san));
fprintf(fid,'FI:  %f\n',FI);
fprintf(fid,'SI:  %f\n',SI);
fprintf(fid,'GRI: %f\n',GRI);
fprintf(fid,'noisy psnr = %f\n',psnrNoisy);
fprintf(fid,'filtered psnr = %f\n',psnrFiltered);
fprintf(fid,'psnr improvement = %f\n',psnrFiltered-psnrNoisy);
fprintf(fid,'noisy img/sig fid met = %f\n',isfmNoisy);
fprintf(fid,'filtered img/sig fid met = %f\n',isfmFiltered);
fprintf(fid,'img/sig fid met improvement = %f\n\n',isfmFiltered-isfmNoisy);
fclose(fid);

%%% Plot results.
start=0; width=260; bottom=0; height=40; 
figure, plot(simArr); axis([start width bottom height]); title('');
figure, plot(simArrNoisy); axis([start width bottom height]); title('');
%figure, plot(filtered); axis([start width bottom height]); title('filtered');
figure, plot(filtered); axis([start width bottom height]); title(filterNameStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Inorm=normalise(I,maxVal);
mx=max(I(:)); mn=min(I(:));
Inorm=(maxVal*(I-mn))/mx-mn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FI=fidelity(org,nos);
% Calculating the fidelity index.
FI=exp(-1*((sum(sum(abs(org-nos))))/(prod(size(org)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yRugStore=dAlphaY(original,recovered);
% Calculating the rugosity in the y-direction.
yRugStore=yDirectionRugosity(original,recovered);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yRugStore=yDirectionRugosity(originalIn,recoveredIn);
% original=clean signal, recovered=denoised(original+gaussian noise).
% This function operates upon the rows of an image there-by finding the
% normal vectors in the y-direction. Each normal vector in the y-direction
% is stored as a row vector, the normal vector for the first row stored in 
% yRugStore(1,:) and the Xth normal vector stored in yRugStore(X,:);
diff_angle_original=[];
diff_angle_recovered=[];
[no_rows,no_cols]=size(originalIn);
yRugStore=zeros(no_rows,no_cols-2);
for i=1:no_rows,
	% Grab rows from the input images.
	original=originalIn(i,:);
	recovered=recoveredIn(i,:);
	for j=1:no_cols-2,
		% Gradient of original signal at point j.
		m_original=original(j)-original(j+1);
		% Find normal vector to original signal at point j.
		n_original=findNormalVector(j,original(j),m_original,no_cols);
		% Gradient of original signal at point j+1.
		m_original_inc=original(j+1)-original(j+2);
		% Find normal vector to original signal at point j+1.
		n_original_inc=findNormalVector(j+1,original(j+1),m_original_inc,no_cols);
		% Calculate the angle between the normal vectors passing through the original signal 
		% at points j and j+1.
		diff_angle_original(1,j)=angleBetween2Vectors(n_original,n_original_inc);
		% Gradient of recovered signal at point j.
		m_recovered=recovered(j)-recovered(j+1);
		% Find normal vector to recovered signal at point j.
		n_recovered=findNormalVector(j,recovered(j),m_recovered,no_cols);
		% Gradient of recovered signal at point j+1.
		m_recovered_inc=recovered(j+1)-recovered(j+2);
		% Find normal vector to recovered signal at point j+1.
		n_recovered_inc=findNormalVector(j+1,recovered(j+1),m_recovered_inc,no_cols);
		% Calculate the angle between the normal vectors passing through the recovered signal 
		% at points j and j+1.
		diff_angle_recovered(1,j)=angleBetween2Vectors(n_recovered,n_recovered_inc);
	end; clear j;
	% Creating a mask indicating the locations where the angle between adjacent normals in the 
	% recovered signal is >= the angle between adjacent normals in the original signal.
	mask=abs(diff_angle_recovered) >= abs(diff_angle_original);
	% Calculate the vertical rugosity in the y-direction.
	rugosity_Ydirection=(diff_angle_recovered-diff_angle_original).*mask;
	% Save this rugosity (in y-direction) vector.
	yRugStore(i,:)=rugosity_Ydirection;
	% Reset variables to be reused.
	diff_angle_original=[];
	diff_angle_recovered=[];
end; clear i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=findNormalVector(x1,y1,m_tan,no_cols);
% Finding the normal vector passing through point (x1,y1) with
% a gradient equal to the negative reciprocal of the gradient of 
% the tangent (m_tan). The normal vector is a straight line of form
% y-y1=m(x-x1).
X=1:no_cols-2;
if m_tan == 0
	m_tan=1e-5;
end
m=-1/m_tan;
y=(m*X)-(m*x1)+y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radians=angleBetween2Vectors(v1,v2);
% Using the dot product to find the angle between two vectors.
dotProdV1V2=sum(v1.*v2);
normV1=sqrt(sum(v1.*v1)); 
normV2=sqrt(sum(v2.*v2)); 
radians=acos(dotProdV1V2/(normV1*normV2));
degrees=(radians*180)/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PSNR = psnr(corrupted, original)
[no_rows,no_cols]=size(original);
MSE = sum(((original-corrupted).^2)) / (no_rows*no_cols);
PSNR = 10*log10(255^2/MSE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isfm=imgSigFidMet(corrupted,original);
[no_rows,no_cols]=size(original);
top=((original-corrupted).^2);
bottom=(original.^2);
isfm=1-(sum(top(:))/sum(bottom(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iwt=wavfilt_ST_PO_log(noisy,thr);
noisy=log(noisy+1);
qmf=MakeONFilter('Daubechies',4);
cl=4;
fwt=FWT_PO(noisy,cl,qmf);
fwt=SoftThresh(fwt,thr);
iwt=IWT_PO(fwt,cl,qmf);
iwt=(exp(iwt))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iwt=wavfilt_HT_PO_log(noisy,thr);
noisy=log(noisy+1);
qmf=MakeONFilter('Daubechies',4);
cl=4;
fwt=FWT_PO(noisy,cl,qmf);
fwt=HardThresh(fwt,thr);
iwt=IWT_PO(fwt,cl,qmf);
iwt=(exp(iwt))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inv=wavfilt_atrousSTlog(noisy,thr);
noisy=log(noisy+1);
simGAUS=zeros(size(noisy));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=FAT_3M(simGAUS,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=simGAUS-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
std0=std(w0);
std1=std(w1);
std2=std(w2);
std3=std(w3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=FAT_3M(noisy,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=noisy-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
tw0=SoftThresh(w0,thr);
tw1=SoftThresh(w1,thr*rat1);
tw2=SoftThresh(w2,thr*rat2);
tw3=SoftThresh(w3,thr*rat3);
inv=s3+tw3+tw2+tw1;
inv=(exp(inv))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inv=wavfilt_atrousHTlog(noisy,thr);
noisy=log(noisy+1);
simGAUS=zeros(size(noisy));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=FAT_3M(simGAUS,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=simGAUS-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
std0=std(w0);
std1=std(w1);
std2=std(w2);
std3=std(w3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=FAT_3M(noisy,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=noisy-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
tw0=HardThresh(w0,thr);
tw1=HardThresh(w1,thr*rat1);
tw2=HardThresh(w2,thr*rat2);
tw3=HardThresh(w3,thr*rat3);
inv=s3+tw3+tw2+tw1+tw0;
inv=(exp(inv))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inv=wavfilt_atrousSSHHTlog(noisy,thr);
noisy=log(noisy+1);
simGAUS=zeros(size(noisy));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=FAT_3M(simGAUS,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=simGAUS-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
std0=std(w0);
std1=std(w1);
std2=std(w2);
std3=std(w3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=FAT_3M(noisy,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=noisy-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
tw0=SoftThresh(w0,thr);
tw1=SoftThresh(w1,thr*rat1);
tw2=HardThresh(w2,thr*rat2);
tw3=HardThresh(w3,thr*rat3);
inv=s3+tw3+tw2+tw1+tw0;
inv=(exp(inv))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inv=wavfilt_atrousSHHHTlog(noisy,thr);
noisy=log(noisy+1);
simGAUS=zeros(size(noisy));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=FAT_3M(simGAUS,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=simGAUS-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
std0=std(w0);
std1=std(w1);
std2=std(w2);
std3=std(w3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=FAT_3M(noisy,0);
s1=FAT_3M(s0,1);
s2=FAT_3M(s1,2);
s3=FAT_3M(s2,3);
w0=noisy-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;
tw0=SoftThresh(w0,thr);
tw1=HardThresh(w1,thr*rat1);
tw2=HardThresh(w2,thr*rat2);
tw3=HardThresh(w3,thr*rat3);
inv=s3+tw3+tw2+tw1+tw0;
inv=(exp(inv))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averaged=cycleSpin(row,th);
qmf=MakeONFilter('Daubechies',4);
rowStore=zeros(size(row));
fprintf('cycle spinning ...\n');

for i=1:length(row),
	rowP=shiftForwardByX(row,i-1);
	%rowP=filterSoft(rowP,qmf,th,'th0+'); % 'th0-' || 'th0+'	
	rowP=filterHard(rowP,qmf,th,'th0+'); % 'th0-' || 'th0+'
	rowP=shiftBackByX(rowP,i-1);
	rowStore(i,:)=rowP;
end; clear i;

averaged=sum(rowStore) / length(row); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filteredRow=filterSoft(row,qmf,th,th0Str);
row=log(row+1);
[l0,h0]=WAVDOWN_1D(row,qmf);
[l1,h1]=WAVDOWN_1D(l0,qmf);
[l2,h2]=WAVDOWN_1D(l1,qmf);
[l3,h3]=WAVDOWN_1D(l2,qmf);

clear l0 l1 l2;

h3=SoftThresh(h3,th);
h2=SoftThresh(h2,th);
h1=SoftThresh(h1,th);
h0=SoftThresh(h0,th);

if strcmp(th0Str,'th0-')==1
	h0=zeros(size(h0));
elseif strcmp(th0Str,'th0+')==1
	h0=h0;
end

l2r=WAVUP_1D(l3,h3,qmf);
l1r=WAVUP_1D(l2r,h2,qmf);
l0r=WAVUP_1D(l1r,h1,qmf);
filteredRow=WAVUP_1D(l0r,h0,qmf);
filteredRow=(exp(filteredRow))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filteredRow=filterHard(row,qmf,th,th0Str);
row=log(row+1);
[l0,h0]=WAVDOWN_1D(row,qmf);
[l1,h1]=WAVDOWN_1D(l0,qmf);
[l2,h2]=WAVDOWN_1D(l1,qmf);
[l3,h3]=WAVDOWN_1D(l2,qmf);

clear l0 l1 l2;

h3=HardThresh(h3,th);
h2=HardThresh(h2,th);
h1=HardThresh(h1,th);
h0=HardThresh(h0,th);

if strcmp(th0Str,'th0-')==1
	h0=zeros(size(h0));
elseif strcmp(th0Str,'th0+')==1
	h0=h0;
end

l2r=WAVUP_1D(l3,h3,qmf);
l1r=WAVUP_1D(l2r,h2,qmf);
l0r=WAVUP_1D(l1r,h1,qmf);
filteredRow=WAVUP_1D(l0r,h0,qmf);
filteredRow=(exp(filteredRow))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new=shiftForwardByX(row,X);
new=zeros(size(row));
new(X+1:length(row))=row(1:length(row)-X);
new(1:X)=row(length(row)-X+1:length(row));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new=shiftBackByX(row,X);
new=zeros(size(row));
new=shiftForwardByX(row,length(row)-X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanOut=datasieve_STlog(noisyIn,th);
noisyIn=log(noisyIn+1);
simGAUS=zeros(size(noisyIn));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=MED_M(simGAUS,3);
s1=MED_M(s0,5);
s2=MED_M(s1,7);
s3=MED_M(s2,9);
d0=simGAUS-s0;
d1=s0-s1;
d2=s1-s2;
d3=s2-s3;
std0=std(d0);
std1=std(d1);
std2=std(d2);
std3=std(d3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=MED_M(noisyIn,3);
s1=MED_M(s0,5);
s2=MED_M(s1,7);
s3=MED_M(s2,9);
d0=noisyIn-s0;
d1=s0-s1;
d2=s1-s2;
d3=s2-s3;
td0=SoftThresh(d0,th);
td1=SoftThresh(d1,th*rat1);
td2=SoftThresh(d2,th*rat2);
td3=SoftThresh(d3,th*rat3);
cleanOut=s3+td3+td2+td1+td0;
cleanOut=(exp(cleanOut))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanOut=datasieve_HTlog(noisyIn,th);
noisyIn=log(noisyIn+1);
simGAUS=zeros(size(noisyIn));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=MED_M(simGAUS,3);
s1=MED_M(s0,5);
s2=MED_M(s1,9);
s3=MED_M(s2,17);
d0=simGAUS-s0;
d1=s0-s1;
d2=s1-s2;
d3=s2-s3;
std0=std(d0);
std1=std(d1);
std2=std(d2);
std3=std(d3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=MED_M(noisyIn,3);
s1=MED_M(s0,5);
s2=MED_M(s1,9);
s3=MED_M(s2,17);
d0=noisyIn-s0;
d1=s0-s1;
d2=s1-s2;
d3=s2-s3;
td0=HardThresh(d0,th);
td1=HardThresh(d1,th*rat1);
td2=HardThresh(d2,th*rat2);
td3=HardThresh(d3,th*rat3);
cleanOut=s3+td3+td2+td1+td0;
cleanOut=(exp(cleanOut))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanOut=datasieve_SSHHTlog(noisyIn,th);
noisyIn=log(noisyIn+1);
simGAUS=zeros(size(noisyIn));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	% From 'type imnoise'.
s0=MED_M(simGAUS,3);
s1=MED_M(s0,5);
s2=MED_M(s1,9);
s3=MED_M(s2,17);
d0=simGAUS-s0;
d1=s0-s1;
d2=s1-s2;
d3=s2-s3;
std0=std(d0);
std1=std(d1);
std2=std(d2);
std3=std(d3);
rat1=std1/std0;
rat2=std2/std0;
rat3=std3/std0;
s0=MED_M(noisyIn,3);
s1=MED_M(s0,5);
s2=MED_M(s1,9);
s3=MED_M(s2,17);
d0=noisyIn-s0;
d1=s0-s1;
d2=s1-s2;
d3=s2-s3;
td0=SoftThresh(d0,th);
td1=SoftThresh(d1,th*rat1);
td2=HardThresh(d2,th*rat2);
td3=HardThresh(d3,th*rat3);
cleanOut=s3+td3+td2+td1+td0;
cleanOut=(exp(cleanOut))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanOut=datasieve_SHHHTlog(noisyIn,th);
noisyIn=log(noisyIn+1);
simGAUS=zeros(size(noisyIn));
simGAUS=simGAUS+sqrt(1)*randn(size(simGAUS))+0;	
