function varargout = first_projet(varargin)
% FIRST_PROJET MATLAB code for first_projet.fig
%      FIRST_PROJET, by itself, creates a new FIRST_PROJET or raises the existing
%      singleton*.
%
%      H = FIRST_PROJET returns the handle to a new FIRST_PROJET or the handle to
%      the existing singleton*.
%
%      FIRST_PROJET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIRST_PROJET.M with the given input arguments.
%
%      FIRST_PROJET('Property','Value',...) creates a new FIRST_PROJET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before first_projet_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to first_projet_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help first_projet

% Last Modified by GUIDE v2.5 11-Feb-2024 22:46:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @first_projet_OpeningFcn, ...
                   'gui_OutputFcn',  @first_projet_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before first_projet is made visible.
function first_projet_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to first_projet (see VARARGIN)

% Choose default command line output for first_projet
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes first_projet wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = first_projet_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Fichier_Callback(hObject, eventdata, handles)
% hObject    handle to Fichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Bruit_Callback(hObject, eventdata, handles)
% hObject    handle to Bruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Transformation_Callback(hObject, eventdata, handles)
% hObject    handle to Transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Filtrage_Callback(hObject, eventdata, handles)
% hObject    handle to Filtrage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to Histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
d = length(size(image));
if d==3
    I = rgb2gray(image);
elseif d==2
    I = image;
end
axes(handles.imgO);
imshow(I);
[nl,nc] = size(I);
v = double(I);
vec = 1:256;
l =0;
for k=0:255 
    for i=1:nl
        for j=1:nc
            if v(i,j)==k 
               l=l+1;
            end
        end
    end
    vec(k+1)=l;
    l=0;
end
axes(handles.imgT);
plot(vec);


% --------------------------------------------------------------------
function Inverser_Callback(hObject, eventdata, handles)
% hObject    handle to Inverser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
image = double(image);
[l,c] = size(image);
image = double(image);
v = image;
for i=1:l
    for j=1:c
        v(i,j) = -double(image(i,j))+255;
    end
end
v = uint8(v);
axes(handles.imgT);
imshow(v);


% --------------------------------------------------------------------
function Lumineuse_Callback(hObject, eventdata, handles)
% hObject    handle to Lumineuse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[l,c] = size(image);
v = image;
for i=1:l
    for j=1:c
        pix = image(i,j)+50;
        if (pix>255)
            pix=255;
        else
            if(pix<0)
                pix = 0;
            end
        end
        v(i,j) = pix;
    end
end
v = uint8(v);
axes(handles.imgT);
imshow(v);


% --------------------------------------------------------------------
function Contrast_Callback(hObject, eventdata, handles)
% hObject    handle to Contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
image = double(image);
[l,c] = size(image);
v = image;
for i=1:l
    for j=1:c
        fpixel = (image(i,j)-128)*5 +128;
        if (fpixel>255)
            fpixel=255;
        else
            if (fpixel<0)
                fpixel = 0;
            end
        end
        v(i,j) = fpixel;
    end
end
v = uint8(v);
axes(handles.imgT);
imshow(v);


% --------------------------------------------------------------------
function Gaussian_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'gaussian', 0,0.01);
%Imgb=handles.Imgb;
axes(handles.imgT)
imshow(Imgb);
handles.courant_data=Imgb;
%subImage(handles.Imgb);
handles.output=hObject;
guidata(hObject,handles);


% --------------------------------------------------------------------
function Poivre_sel_Callback(hObject, eventdata, handles)
% hObject    handle to Poivre_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'Salt & Pepper', 0.01);
%Imgb=handles.Imgb;
axes(handles.imgT)
%
imshow(Imgb);
handles.courant_data=Imgb;
%subImage(handles.Imgb);
handles.output=hObject;
guidata(hObject,handles);


% --------------------------------------------------------------------
function Poissan_Callback(hObject, eventdata, handles)
% hObject    handle to Poissan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Ouvrir_Callback(hObject, eventdata, handles)
% hObject    handle to Ouvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.imgO)
handles.courant_data = handles.ima;
imshow(handles.courant_data);
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Enregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to Enregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Enregistrer votre image ...');
imwrite(image, sprintf('%s',path,file),'png');


% --------------------------------------------------------------------
function Quitter_Callback(hObject, eventdata, handles)
% hObject    handle to Quitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --------------------------------------------------------------------
function Gauss_3_3_Callback(hObject, eventdata, handles)
% hObject    handle to Gauss_3_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/16)*[1 2 1 ;2 4 2 ; 1 2 1];
for x = 2 : n-1
    for y = 2 : m-1
    f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Gauss_5_5_Callback(hObject, eventdata, handles)
% hObject    handle to Gauss_5_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/256)*[1 4 6 4 1 ; 4 16 24 16 4 ; 6 24 36 24 6 ; 4 16 24 16 4 ; 1 4 6 4 1];
for x = 3 : n-2
    for y = 3 : m-2
  f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Laplacien_Callback(hObject, eventdata, handles)
% hObject    handle to Laplacien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%image=imnoise(imageO,'salt & pepper', 0.05);
image = double(image);
%b=image;
[n,m]=size(image);
b=zeros(n,m);
%M1=[0 1 0;1 -4 1;0 1 0];
M1=[-1 -1 -1;-1 8 -1;-1 -1 -1];
for i=2:n-1
    for j=2:m-1
        V=image(i-1:i+1,j-1:j+1);
        S=V.*M1;
        b(i,j)=sum(S(:));
    end
end
b=uint8(b);
axes(handles.imgT);
imshow(b);


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
%image=rgb2gray(image);

output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
%maskhor = [-1,0,1;-b,0,b;-1,0,1]; 
%maskver = [-1,-b,-1;0,0,0;1,b,1];
%b=2 pour sobel
maskhor = [-1,0,1;-2,0,2;-1,0,1]; 
maskver = [-1,-2,-1;0,0,0;1,2,1];

for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
   for j=4:(n-3) 
output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 
axes(handles.imgT);
imshow(output);



% --------------------------------------------------------------------
function Robert_Callback(hObject, eventdata, handles)
% hObject    handle to Robert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Conique_Callback(hObject, eventdata, handles)
% hObject    handle to Conique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[ n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[0 0 1 0 0 ; 0 2 2 2 0 ; 1 2 5 2 1 ; 0 2 2 2 0 ; 0 0 1 0 0];
for x = 3 : n-2
    for y = 3 : m-2
       f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);



% --------------------------------------------------------------------
function Median_Callback(hObject, eventdata, handles)
% hObject    handle to Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image=double(image);
[n,m]=size(image);
img=image;
for i=2:n-1
    for j=2:m-1
       fenetre=image(i-1:i+1,j-1:j+1);
       v=[fenetre(1,:) fenetre(2,:) fenetre(3,:)];
       sort(v);
       a=median(v);
       img(i,j)=a;
    end
end
b=uint8(img);
handles.ima_traite = b;
axes(handles.imgT);
imshow(b);
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Moy_3_3_Callback(hObject, eventdata, handles)
% hObject    handle to Moy_3_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/9)*[1 1 1 ;1 1 1  ; 1 1 1];
for x = 2 : n-1
    for y = 2 : m-1
    f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Moy_5_5_Callback(hObject, eventdata, handles)
% hObject    handle to Moy_5_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ];
for x = 3 : n-2
    for y = 3 : m-2
  f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Filtre_frequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to Filtre_frequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Niveau_de_gris_Callback(hObject, eventdata, handles)
% hObject    handle to Niveau_de_gris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ima = handles.courant_data;
d = length(size(ima));
if d==3
    imagray=rgb2gray(ima);
elseif d==2
    imagray=ima;
end
axes(handles.imgT);
imshow(imagray);


% --------------------------------------------------------------------
function Binairisation_Callback(hObject, eventdata, handles)
% hObject    handle to Binairisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I4 = handles.courant_data;
m0=1;
m1=mean2(I4); m2 = mean2(I4.^2); m3=mean2(I4.^3);
C1=(m3-(m1*m2))/(m2-m1);
C0=(-m2-(C1*m1))/m0;
z1 = (-C1-sqrt(C1^2-4*C0))/2;
z2 = (-C1+sqrt(C1^2-4*C0))/2;
seuil = (z1+z2)/2;
bin = (I4>seuil)*255;
text(2,10,num2str(seuil));
axes(handles.imgO);
imshow(I4);
axes(handles.imgT);
handles.ima_traite = bin;
imshow(handles.ima_traite);
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Fpb_Callback(hObject, eventdata, handles)
% hObject    handle to Fpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data;
 F=fftshift(fft2(I)); 
% %calcul de la taille de l'image; 
 M=size(F,1); 
 N=size(F,2); 
 P=size(F,3); 
 H0=zeros(M,N); 
 D0=1; 
 M2=round(M/2); 
 N2=round(N/2); 
 H0(M2-D0:M2+D0,N2-D0:N2+D0)=1; 
 for i=1:M 
  for j=1:N 
      G(i,j)=F(i,j)*H0(i,j); 
 end 
 end 
 g=ifft2(G); 
 imshow(abs(g),[0,255]);



% --------------------------------------------------------------------
function Fpbb_Callback(hObject, eventdata, handles)
% hObject    handle to Fpbb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data;
%I = imread('eight.tif');

F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3);

H0=zeros(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H0(M2-D0:M2+D0,N2-D0:N2+D0)=1; 

n=3; 

for i=1:M 
for j=1:N 
%H(i,j)=1/(1+(H0(i,j)/D0)^(2*n)); 
G(i,j)=F(i,j)*H0(i,j); 
end 
end 

g=ifft2(G); 

%subplot(1,2,1);imshow(I);title('image originale'); 
%subplot(1,2,2);
imshow(abs(g),[0,255]);%title('image filtrée'); 


% --------------------------------------------------------------------
function Fph_Callback(hObject, eventdata, handles)
% hObject    handle to Fph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
%charge; 
F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 

H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 
for i=1:M 
for j=1:N 
G(i,j)=F(i,j)*H1(i,j); 
end 
end 
g=ifft2(G); 
%subplot(1,2,1);imshow(I);title('image originale'); 
%subplot(1,2,2);
imshow(255-abs(g),[0,255]);



% --------------------------------------------------------------------
function Fphb_Callback(hObject, eventdata, handles)
% hObject    handle to Fphb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;

F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 

H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 

n=3; 

for i=1:M 
for j=1:N 
H(i,j)=1/(1+(H1(i,j)/D0)^(2*n)); 
G(i,j)=F(i,j)*H(i,j); 
end 
end 

g=ifft2(G); 

%subplot(1,2,1);imshow(I);title('image originale'); 
%subplot(1,2,2);
imshow(255-abs(g),[0,255]);%title('image filtrée');



% --------------------------------------------------------------------
function Gradient_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image = double(image);
%image=rgb2gray(image);

[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [0,0,0;-1,0,1;0,0,0]; 
maskver = [0,-1,0;0,0,0;0,1,0];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3         
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);            
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
for j=4:(n-3)       
    output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j));
end 
end 
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 

%b=uint8(b);
axes(handles.imgT);
imshow(output);


% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
%image=rgb2gray(image);

output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
%maskhor = [-1,0,1;-b,0,b;-1,0,1]; 
%maskver = [-1,-b,-1;0,0,0;1,b,1];
%b=1 pour prewitt
maskhor = [-1,0,1;-1,0,1;-1,0,1]; 
maskver = [-1,-1,-1;0,0,0;1,1,1];

for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
   for j=4:(n-3) 
output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 
axes(handles.imgT);
imshow(output);



% --------------------------------------------------------------------
function Pyramidal_Callback(hObject, eventdata, handles)
% hObject    handle to Pyramidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/81)*[1 2 3 2 1 ; 2 4 6 4 2 ; 3 6 9 6 3 ; 2 4 6 4 2 ; 1 2 3 2 1];
for x = 3 : n-2
    for y = 3 : m-2
          f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function morphologie_Callback(hObject, eventdata, handles)
% hObject    handle to morphologie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function point_dinterer_Callback(hObject, eventdata, handles)
% hObject    handle to point_dinterer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function hough_Callback(hObject, eventdata, handles)
% hObject    handle to hough (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
siz = size(image);
rl = ceil(sqrt(siz(1)^2+siz(2)^2));

h = zeros(rl,360);
theta = 360;

% accumulateur
for x = 1:siz(1)
    for y = 1:siz(2)
        if image(x,y)==1
            for theta = 1:360
                r = round(x*cos(theta*pi/180) + y*sin(theta*pi/180));
                if r<1
                    r = 1;
                    h(r,theta) = h(r,theta);
                else
                    h(r,theta) = h(r,theta) + 1;
                end
            end
        end
    end
end
subplot(1,2,2)
imshow(h)




% --------------------------------------------------------------------
function susan_Callback(hObject, eventdata, handles)
% hObject    handle to susan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
im=handles.courant_data;
d = length(size(im));
% im = imnoise(im,'speckle',0.02);
if d==3
    image=double(rgb2gray(im));
elseif d==2
    image=double(im);
end
subplot(1,2,1)

imshow(im)
[n,m]=size(image);
% =============================données=====================================
rayon=2;
alpha=50;
r=2;
alpha=alpha/100;
% ========================génerateur de mask===============================
mask=zeros(2*rayon+1);
b=ones(rayon+1);
for i=1:rayon+1
    for j=1:rayon+1
        if (rayon==1)
           if(j>i)
            b(i,j)=0;
           end
         else
             if(j>i+1)
            b(i,j)=0;
         end
        end
    end
end
mask(1:rayon+1,rayon+1:2*rayon+1)=b;
mask(1:rayon+1,1:rayon+1)=rot90(b);
mask0=mask;
mask0=flipdim(mask0,1);
mask=mask0+mask;
mask(rayon+1,:)=mask(rayon+1,:)-1;
% ==========================réponse maximale===============================
max_reponse=sum(sum(mask));
% =====================balayage de toute l'image===========================
f=zeros(n,m);
for i=(rayon+1):n-rayon
    for j=(rayon+1):m-rayon
  
          image_courant=image(i-rayon:i+rayon,j-rayon:j+rayon);

    image_courant_mask=image_courant.*mask;

         inteniste_cental= image_courant_mask(rayon+1,rayon+1);
         s=exp(-1*(((image_courant_mask-inteniste_cental)/max_reponse).^6));
       somme=sum(sum(s));
%   si le centre du mask est un 0 il faut soustraire les zeros des filtres
                if (inteniste_cental==0)
                    somme=somme-length((find(mask==0)));
                end       
         f(i,j)=somme;           
     end
end
% =============selection et seuillage des points d'interét=================
ff=f(rayon+1:n-(rayon+1),rayon+1:m-(rayon+1));
minf=min(min(ff));
maxf=max(max(f));
fff=f;
d=2*r+1;
temp1=round(n/d);
if (temp1-n/d)<0.5 &(temp1-n/d)>0
    temp1=temp1-1;
end
temp2=round(m/d);
if (temp2-m/d)<0.5 &(temp2-m/d)>0
    temp2=temp2-1;
end
fff(n:temp1*d+d,m:temp2*d+d)=0;
temp1
temp2

fff;

sizefff=size(fff)
for i=(r+1):d:temp1*d+d
for j=(r+1):d:temp2*d+d
     window=fff(i-r:i+r,j-r:j+r);
     window0=window;
   [xx,yy]=find(window0==0);
               for k=1:length(xx)
                 window0(xx(k),yy(k))=max(max(window0));
               end
               minwindow=min(min(window0));
[y,x]=find(minwindow~=window & window<=minf+alpha*(maxf-minf) & window>0);
[u,v]=find(minwindow==window);
if length(u)>1
             for l=2:length(u) 
                     fff(i-r-1+u(l),j-r-1+v(l))=0   ;       
             end
end
if length(x)~=0
        for l=1:length(y) 
            fff(i-r-1+y(l),j-r-1+x(l))=0   ;
        end
end
fff;
end
end


seuil=minf+alpha*(maxf-minf);

    [u,v]=find(minf<=fff & fff<=seuil );
% ==============affichage des resultats====================================
subplot(1,2,2)
imshow(im)
hold on
plot(v,u,'.r','MarkerSize',10)
nombre_de_point_dinteret=length(v)


% --------------------------------------------------------------------
function harris_Callback(hObject, eventdata, handles)
% hObject    handle to harris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img=handles.courant_data;
%==========================================================================
if(size(img,3)==3)
    display('l''image est en couleur')
    img=rgb2gray(img);
end
%==========================================================================

lambda=0.04;
sigma=1; seuil=200; r=6; w=5*sigma;
[m,n]=size(img)
imd=double(img);
dx=[-1 0 1
    -2 0 2
    -1 0 1]; % derivée horizontale : filtre de Sobel
dy=dx'; % derivée verticale : filtre de Sobel
% dx=[-1 1];
% dy=dx';
g = fspecial('gaussian',max(1,fix(w)), sigma);
Ix=conv2(imd,dx,'same');
Iy=conv2(imd,dy,'same');
Ix2=conv2(Ix.^2, g, 'same');
Iy2=conv2(Iy.^2, g, 'same');
Ixy=conv2(Ix.*Iy, g,'same');
% M=[Ix2 Ixy
%     Ixy Iy2];
detM=Ix2.*Iy2-Ixy.^2;
trM=Ix2+Iy2;
R=detM-lambda*trM.^2;
%==========================================================================
R1=(1000/(1+max(max(R))))*R;
%==========================================================================          
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
    R1(u(k),v(k))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);
for i=r+1:m1-r
    for j=r+1:n1-r
        fenetre=R11(i-r:i+r,j-r:j+r);
        ma=max(max(fenetre));
        if fenetre(r+1,r+1)<ma
            R11(i,j)=0;
        end
    end
end
axes(handles.imgT);
imshow(img);
handles.ima_traite = img;
handles.output = hObject;
guidata(hObject, handles);
hold on
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
nb=length(x)
plot(y,x,'.r')
title('detection des points d''interets')


% --------------------------------------------------------------------
function erosion_Callback(hObject, eventdata, handles)
% hObject    handle to erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger une image en niveaux de gris
image = handles.courant_data;

% Convertir l'image en niveaux de gris si elle est en couleur
if size(image, 3) == 3
    image = rgb2gray(image);
end

% Définir l'élément structurant pour l'érosion (dans cet exemple, un carré de taille 3x3)
se = strel('square', 3);

% Effectuer l'érosion de l'image
image_eroded = imerode(image, se);

% Afficher l'image originale et l'image érodée
subplot(1, 2, 2);
imshow(image_eroded);
title('Image érodée');


% --------------------------------------------------------------------
function dilatation_Callback(hObject, eventdata, handles)
% hObject    handle to dilatation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger une image en niveaux de gris
image = handles.courant_data;
% Convertir l'image en niveaux de gris si elle est en couleur
if size(image, 3) == 3
    image = rgb2gray(image);
end

% Définir l'élément structurant pour la dilatation (dans cet exemple, un carré de taille 3x3)
se = strel('square', 3);

% Effectuer la dilatation de l'image
image_dilated = imdilate(image, se);

% Afficher l'image originale et l'image dilatée
subplot(1, 2, 2);
imshow(image_dilated);
title('Image dilatée');


% --------------------------------------------------------------------
function ouverture_Callback(hObject, eventdata, handles)
% hObject    handle to ouverture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger une image en niveaux de gris
image = handles.courant_data;
% Convertir l'image en niveaux de gris si elle est en couleur
if size(image, 3) == 3
    image = rgb2gray(image);
end

% Définir l'élément structurant pour l'ouverture (dans cet exemple, un carré de taille 3x3)
se = strel('square', 3);

% Effectuer l'ouverture de l'image
image_opened = imopen(image, se);

% Afficher l'image originale et l'image ouverte
subplot(1, 2, 2);
imshow(image_opened);
title('Image ouverte');


% --------------------------------------------------------------------
function fermeture_Callback(hObject, eventdata, handles)
% hObject    handle to fermeture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger une image en niveaux de gris
image =handles.courant_data;

% Convertir l'image en niveaux de gris si elle est en couleur
if size(image, 3) == 3
    image = rgb2gray(image);
end

% Définir l'élément structurant pour la fermeture (dans cet exemple, un carré de taille 3x3)
se = strel('square', 3);

% Effectuer la fermeture de l'image
image_closed = imclose(image, se);

% Afficher l'image originale et l'image fermée
subplot(1, 2, 2);
imshow(image_closed);
title('Image fermée');
