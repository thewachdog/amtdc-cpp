clear all
% The program is generated for the dynamic interaction of pantograph and
% OHE system. the length of the element will be constant for this as 0.25
% meters. The code is generalized for any number of droppers per span. The
% input required from users are: dropper locations and the maximum presag
% in the contact wire
%% System parameters
% rho=linear density C: contact wire, M: Messenger wire, D: droppers.
% A= cross-sectional area  C: contact wire, M: Messenger wire, D: droppers.
% E=Young's modulus, C: contact wire, M: Messenger wire, D: droppers.
% I=Area moment of Inertia, C: contact wire, M: Messenger wire, D: droppers.
% T=Tension in the wire, C: contact wire, M: Messenger wire

% Contact wire
rho_C =1.35;A_C = 150*10^-6;E_C = 100*10^9;I_C = 1.7905*10^-9; EA_C = E_C*A_C; EI_C = E_C*I_C; T_C=-22000;

% Mesenger wire
rho_M=1.080; A_M=120*10^-6; E_M=97*10^9; I_M=1.1459*10^-9; EA_M=E_M*A_M; EI_M=E_M*I_M; T_M=-16000;

% Droppers
rho_D=0.117; EA_D=200*10^3; EI_D=0.46;

No_span=2;
% No_span=input('Number of spans:  '); % Please enter the number of spans you want to simulate

% contact wire support
EA_S_CW=14167000; EI_S_CW=0; 

% Messenger wire support
EA_S_MW= 500000; EI_S_MW=0;

lengthSpan=55; % Length of one span
T_c=-T_C;
T_m=-T_M;
w_c=rho_C*9.8; % weight density of the contact wire
w_m=rho_M*9.8; % weight density of the messenger wire
w_dm=0.195*9.8; %weight of dropper clamp on messenger wire
w_dc=0.165*9.8; % weight of dropper clamp on contact wire

No_droppers=9;
% No_droppers=input('Number of droppers per span: '); % Number of droppers per span

%for i=1:No_droppers
%    sprintf('Enter the location of dropper %d ' ,i)
%    dropper_location(i)=input(':');      % x coordinate of the droppers location 
%end

dropper_location=[4.5 10.25 16 21.75 27.5 33.25 39 44.75 50.5];

% d_c=input("Enter the required presag of contact wire:   "); % Maximum presag on the contact wire
d_c = 0.055;
encumbrance=1.2;

%% Dropper force and length

% presag calculated at each droppper location
D_dropper=(4*d_c/((dropper_location(No_droppers)-dropper_location(1))^2)).*(dropper_location-dropper_location(1)).*((dropper_location(No_droppers)-dropper_location(1))-dropper_location+dropper_location(1));

% Forces at the droppers
F=zeros(No_droppers,1);
f1=dropper_location(1);
f2=(dropper_location(2)-dropper_location(1))/2;
f3=(T_c*(D_dropper(2)-D_dropper(1)))/(w_c*(dropper_location(2)-dropper_location(1)));
F_1_9=(f1+f2+f3)*w_c+w_dc;
F(1)=F_1_9;
F(9)=F_1_9;

for i=2:(No_droppers-1)
    f4=(dropper_location(i+1)-dropper_location(i-1))/2;
    f5=(T_c*(D_dropper(i)-D_dropper(i-1)))/(w_c*(dropper_location(i)-dropper_location(i-1)));
    f6=(T_c*(D_dropper(i+1)-D_dropper(i)))/(w_c*(dropper_location(i+1)-dropper_location(i)));
    F(i)=(f4-f5+f6)*w_c+w_dc;
end

% Reaction force at the messenger wire support
for i=1:No_droppers
    r1(i)=(F(i)*(lengthSpan-dropper_location(i)))/lengthSpan;
end
r2=sum(r1(1:No_droppers));
R_a=((w_m*lengthSpan)/2)+r2;

F_droppers_m=zeros(No_droppers,1);
F_droppers_m(1)=0;
for i=2:No_droppers
    for j=1:i-1
        F_droppers_m1(j)=F(j)*5.75*(i-j);
    end

    F_droppers_m(i)=sum(F_droppers_m1);
    
    if i> ((No_droppers/2)+1)
       F_droppers_m(i)=F_droppers_m((No_droppers+1)-i);
    end
end

% Sag calculation at the messenger wire at the dropper locations
for i=1:No_droppers
    c1=R_a*dropper_location(i);
    c2=(w_m*(dropper_location(i))^2)/2;
    c3=sum(F_droppers_m(i));
    c_m(i)=(c1-c2-c3)/T_m;
        if i> ((No_droppers/2)+1)
       c_m(i)=c_m((No_droppers+1)-i);
    end
end

% Calculation of dropper lengths 
for i=1:No_droppers
    length_dropper(i)=encumbrance-c_m(i)+D_dropper(i);
end
  
%% Node coordinates

EL=0.25; % length of each element in the wire

numberElements=No_span*449+2*No_span+2; % total number of elements (220+220+9) contact wire, messenger wire, droppers

numberNodes_w=221; % number of elements for contact wire and messenger wire for a single span

numberElements_W=(No_span*lengthSpan/EL);% Total Number of elements for contact wire and messneger wire

numberNodes_W=(No_span*lengthSpan/EL)+1;% Total Number of nodes for contact wire and messneger wire


NC_CW=zeros(numberNodes_W,2); % nodal coordinates of the contact wire
NC_MW=zeros(numberNodes_W,2); % % nodal coordinates of the messenger wire

mw_x=0:0.25:55;% Discretization of the wire into finite elements for sigle span
Nodes_effective = (2*numberNodes_W+2*No_span+2);% Total number of nodes of the OHE system

%dropper position and support on the contact wire (x coordinates)
cW_x=zeros((No_droppers+2),1);
for i= 2:(No_droppers+1)
    cW_x(i)=dropper_location(i-1);
    cW_x(No_droppers+2)=lengthSpan;
end
% dropper position and support on the contact wire (y coordinates)
cW_y=zeros(No_droppers+2,1);
for i= 2:(No_droppers+1)
    cW_y(i)=D_dropper(i-1);
end
%Interpolation of the dropper point to find out the configuration of the entire contact wire
CW_Y = interp1(cW_x,cW_y,mw_x);
CW_Y=-CW_Y;

MW_x=cW_x; % messenger wire x coordinates

MW_y=zeros(No_droppers+2,1);% messenger wire y coordinates for dropper positions
for i=1:No_droppers+2
    if i==1||i==No_droppers+2
        MW_y(i)=1.2;
    else
        MW_y(i)=length_dropper(i-1);
    end
end

        
MW_Y = interp1(MW_x,MW_y,mw_x,'spline'); % Interpolation to find out the y coordinates of the messenger wire for finding its initial position


for i=1:numberNodes_w % Assigning the nodal coordinates to the messenger wire
    NC_MW(i,1)=mw_x(i);
end

% Node coordinates for contact wire
for j=1:No_span
    for i=1:221
        if j==1
            NC_CW(i,1)=mw_x(i);
            NC_CW(i,2)=CW_Y(i);
        else
            NC_CW((i+(j-1)*221)-(j-1),1)=mw_x(i)+(j-1)*lengthSpan;
            NC_CW((i+(j-1)*221)-(j-1),2)=CW_Y(i);
        end
    end
end

% Node coordinates for Messenger wire
for j=1:No_span
    for i=1:221
        if j==1
            NC_MW(i,1)=mw_x(i);
            NC_MW(i,2)=MW_Y(i);
        else
            NC_MW((i+(j-1)*221)-(j-1),1)=mw_x(i)+(j-1)*lengthSpan;
            NC_MW((i+(j-1)*221)-(j-1),2)=MW_Y(i);
        end
    end
end

%% Droppers Node coordinates
NC_droppers_single=zeros(2*No_droppers,2);
for i=1:No_droppers
    NC_droppers_single(2*i-1,1)=dropper_location(i);
   NC_droppers_single(2*i,1)=dropper_location(i);
   NC_droppers_single(2*i,2)=length_dropper(i);      
end
%Dropper node coordinates for multiple span
NC_droppers=zeros(No_span*18,2);

for j=1:No_span
    for i=1:18
        if j==1
            NC_droppers(i,1)=NC_droppers_single(i,1);
            NC_droppers(i,2)=NC_droppers_single(i,2);
        else
            
            NC_droppers((i+(j-1)*18),1)=NC_droppers_single(i,1)+(j-1)*55;
            NC_droppers((i+(j-1)*18),2)=NC_droppers_single(i,2);

        end
    end
end
 
%% Supports at the wires
NC_supports_CW = zeros(2*No_span+2,2);
for i=1:2:length(NC_supports_CW)
    NC_supports_CW(i,2)=0;
     NC_supports_CW(i,1)=((i-1)/2)*55;
end
for i=2:2:length(NC_supports_CW)
    NC_supports_CW(i,2)=-1;
     NC_supports_CW(i,1)=((i-2)/2)*55;
end

NC_supports_MW = zeros(2*No_span+2,2);
for i=1:2:length(NC_supports_MW)
    NC_supports_MW(i,2)=1.2;
     NC_supports_MW(i,1)=((i-1)/2)*55;
end
for i=2:2:length(NC_supports_MW)
    NC_supports_MW(i,2)=2.2;
     NC_supports_MW(i,1)=((i-2)/2)*55;
end

%% Combining different matrices
nodeCoordinates=zeros(length(Nodes_effective),2);
% node coordinates for Contact wire
for i=1:length(NC_CW)
    for j=1:2
    nodeCoordinates(i,j)=NC_CW(i,j);
    end
end

% node coordinates for Messenger wire
for i=length(NC_CW)+1:2*length(NC_MW)
    for j=1:2
    nodeCoordinates(i,j)=NC_MW(i-length(NC_CW),j);
    end
end

for i=(2*length(NC_MW)+1):(2*length(NC_MW))+(No_span+1)
        for j=1:2
    nodeCoordinates(i,j)=NC_supports_CW(2*(i-(2*length(NC_MW))),j);
    end
end
 
for i=(2*length(NC_MW))+(No_span+1)+1:(2*length(NC_MW))+(No_span+1)+(No_span+1)
        for j=1:2
    nodeCoordinates(i,j)=NC_supports_MW(2*(i-((2*length(NC_MW))+(No_span+1))),j);
    end
end

%% Connection between different elements

numberNodes_W=(No_span*lengthSpan/EL)+1;
numberElements_W=(No_span*lengthSpan/EL);

% Connection between contact wire elements nodes
for i=1:numberElements_W
    elementNodes(i,1)=i;
    elementNodes(i,2)=i+1;
end

% Connection between messenger wire elements nodes
for i=(numberElements_W+1):(2*numberElements_W)
    elementNodes(i,1)=i+1;
    elementNodes(i,2)=i+2;
end

 % Connection between droppers element nodes
 for i=(2*numberElements_W)+1:(2*numberElements_W)+No_span*9
     elementNodes(i,1)=((NC_droppers(2*(i-(2*numberElements_W))))/EL)+1;
     elementNodes(i,2)=(((NC_droppers(2*(i-(2*numberElements_W))))/EL)+1)+((No_span*lengthSpan)/EL)+1;
 end

%Connections between contact wire support
for i=(No_span*9+2*numberElements_W)+1:(No_span*9+2*numberElements_W)+1+No_span
      elementNodes(i,1)= ((i-((No_span*9+2*numberElements_W)+1))*lengthSpan)/EL+1;
      elementNodes(i,2) =(i-(No_span*9-2));
end
    
% Connection between Messenger wire support element nodes 
for i=(No_span*9+2*numberElements_W)+1+No_span+1:(No_span*9+2*numberElements_W)+2+2*No_span
      elementNodes(i,1)= (((i-((No_span*9+2*numberElements_W)+1+No_span+1))*lengthSpan)/EL)+1+numberNodes_W;
      elementNodes(i,2) =(i-(No_span*9-2));
end

xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

%% Numbner of degrees of freedom and effective number of nodes

GDof = 3*(2*((No_span*lengthSpan/EL)+1)+(2*No_span+2));

Nodes_effective = (2*numberNodes_W+2*No_span+2);

%% external force

% Self weight
force=zeros(GDof,1);
Length_Element=0.25;
% Self weight of the contact wire
for i=(Nodes_effective+1):(Nodes_effective+numberElements_W)
    force(i)=force(i)-(9.81*1.35*Length_Element/4);
end

%% Stiffness formation

stiffness=zeros(GDof);

% Contact wire
for e = 1:numberElements_W
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];
    nn = length(indice);
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    cosa = xa/length_element;
    sena = ya/length_element;
    ll = length_element;
    L = [cosa*eye(2) sena*eye(2) zeros(2); -sena*eye(2) cosa*eye(2) zeros(2); zeros(2,4) eye(2)];
    oneu = [1 -1;-1 1];
    oneu2 = [1 -1;1 -1];
    oneu3 = [1 1;-1 -1];
    oneu4 = [4 2;2 4];

    A=[(4*EI_C/ll)+(4*ll*T_C/30) (2*EI_C/ll)-(ll*T_C/30); (2*EI_C/ll)-(ll*T_C/30)  (4*EI_C/ll)+(4*ll*T_C/30)];
    
    k1 = [EA_C/ll*oneu zeros(2,4);
    zeros(2) ((36*T_C/30*ll)+(12*EI_C/ll^3))*oneu ((3*T_C/30)+(6*EI_C/ll^2))*oneu3;
    zeros(2) ((3*T_C/30)+(6*EI_C/ll^2))*oneu2 A];
    stiffness(elementDof,elementDof) = stiffness(elementDof,elementDof) + L'*k1*L;
end

% Messenger wire
for e = numberElements_W+1:2*numberElements_W
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];
    nn = length(indice);
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    cosa = xa/length_element;
    sena = ya/length_element;
    ll = length_element;
    L = [cosa*eye(2) sena*eye(2) zeros(2); -sena*eye(2) cosa*eye(2) zeros(2); zeros(2,4) eye(2)];
    oneu = [1 -1;-1 1];
    oneu2 = [1 -1;1 -1];
    oneu3 = [1 1;-1 -1];
    oneu4 = [4 2;2 4];
    
    A=[(4*EI_M/ll)+(4*ll*T_M/30) (2*EI_M/ll)-(ll*T_M/30); (2*EI_M/ll)-(ll*T_M/30)  (4*EI_M/ll)+(4*ll*T_M/30)];
    
    k1 = [EA_M/ll*oneu zeros(2,4);
    zeros(2) ((36*T_M/30*ll)+(12*EI_M/ll^3))*oneu ((3*T_M/30)+(6*EI_M/ll^2))*oneu3;
    zeros(2) ((3*T_M/30)+(6*EI_M/ll^2))*oneu2 A];
    stiffness(elementDof,elementDof) = stiffness(elementDof,elementDof) + L'*k1*L;
end

% Dropper
for e = 2*numberElements_W+1:(2*numberElements_W+9*No_span)
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];
    nn = length(indice);
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    cosa = xa/length_element;
    sena = ya/length_element;
    ll = length_element;
    L = [cosa*eye(2) sena*eye(2) zeros(2); -sena*eye(2) cosa*eye(2) zeros(2); zeros(2,4) eye(2)];
    oneu = [1 -1;-1 1];
    oneu2 = [1 -1;1 -1];
    oneu3 = [1 1;-1 -1];
    oneu4 = [4 2;2 4];
         
    T=0;
    
    A=[(4*EI_D/ll)+(4*ll*T/30) (2*EI_D/ll)-(ll*T/30); (2*EI_D/ll)-(ll*T/30)  (4*EI_D/ll)+(4*ll*T/30)];
    
    k1 = [EA_D/ll*oneu zeros(2,4);
        zeros(2) ((36*T/30*ll)+(12*EI_D/ll^3))*oneu ((3*T/30)+(6*EI_D/ll^2))*oneu3; 
        zeros(2) ((3*T/30)+(6*EI_D/ll^2))*oneu2 A];

    stiffness(elementDof,elementDof) = stiffness(elementDof,elementDof) + L'*k1*L;
end

% contact wire support
for e = (2*numberElements_W+9*No_span)+1:(No_span*9+2*numberElements_W)+1+No_span
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];   
    nn = length(indice);
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    cosa = xa/length_element;
    sena = ya/length_element;
    ll = length_element;
    L = [cosa*eye(2) sena*eye(2) zeros(2); -sena*eye(2) cosa*eye(2) zeros(2); zeros(2,4) eye(2)];
    oneu = [1 -1;-1 1];
    oneu2 = [1 -1;1 -1];
    oneu3 = [1 1;-1 -1];
    oneu4 = [4 2;2 4];

    T=0;
    
    A=[(4*EI_S_CW/ll)+(4*ll*T/30) (2*EI_S_CW/ll)-(ll*T/30); (2*EI_S_CW/ll)-(ll*T/30)  (4*EI_S_CW/ll)+(4*ll*T/30)];
    
    k1 = [EA_S_CW/ll*oneu zeros(2,4);
        zeros(2) ((36*T/30*ll)+(12*EI_S_CW/ll^3))*oneu ((3*T/30)+(6*EI_S_CW/ll^2))*oneu3; 
        zeros(2) ((3*T/30)+(6*EI_S_CW/ll^2))*oneu2 A];

    stiffness(elementDof,elementDof) = stiffness(elementDof,elementDof) + L'*k1*L;
end

% Messenger wire support

for e = (No_span*9+2*numberElements_W)+1+No_span+1:(No_span*9+2*numberElements_W)+2+2*No_span
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];    
    nn = length(indice);
    xa = xx(indice(2))-xx(indice(1));
    ya = yy(indice(2))-yy(indice(1));
    length_element = sqrt(xa*xa+ya*ya);
    cosa = xa/length_element;
    sena = ya/length_element;
    ll = length_element;
    L = [cosa*eye(2) sena*eye(2) zeros(2); -sena*eye(2) cosa*eye(2) zeros(2); zeros(2,4) eye(2)];
    oneu = [1 -1;-1 1];
    oneu2 = [1 -1;1 -1];
    oneu3 = [1 1;-1 -1];
    oneu4 = [4 2;2 4];
         
    T=0;
    
    A=[(4*EI_S_MW/ll)+(4*ll*T/30) (2*EI_S_MW/ll)-(ll*T/30); (2*EI_S_MW/ll)-(ll*T/30)  (4*EI_S_MW/ll)+(4*ll*T/30)];
    
    k1 = [EA_S_MW/ll*oneu zeros(2,4);
        zeros(2) ((36*T/30*ll)+(12*EI_S_MW/ll^3))*oneu ((3*T/30)+(6*EI_S_MW/ll^2))*oneu3; 
        zeros(2) ((3*T/30)+(6*EI_S_MW/ll^2))*oneu2 A];

    stiffness(elementDof,elementDof) = stiffness(elementDof,elementDof) + L'*k1*L;
end

%% Mass Matrix formulation

mass = zeros(GDof);

% Contact wire
for e = 1:numberElements_W
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];
xa = xx(indice(2))-xx(indice(1));
ya = yy(indice(2))-yy(indice(1));
length_element = sqrt(xa*xa+ya*ya);
cosa = xa/length_element;
sena = ya/length_element;
ll = length_element;
L = [cosa*eye(2) sena*eye(2) zeros(2);
-sena*eye(2) cosa*eye(2) zeros(2);
zeros(2,4) eye(2)];
oneu = 1/6*[2 1;1 2];
oneu2 = 1/70*[26 9; 9 26];
oneu3 = ll/420*[22 -13; 13 -22];
oneu4 = ll^2/420*[4 -3; -3 4];
m1 = rho_C*ll*[oneu zeros(2,4);
zeros(2) oneu2 oneu3;
zeros(2) oneu3' oneu4];
mass(elementDof,elementDof) = mass(elementDof,elementDof) + L'*m1*L;
end

% Messenger wire
for e = numberElements_W+1:2*numberElements_W
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];
xa = xx(indice(2))-xx(indice(1));
ya = yy(indice(2))-yy(indice(1));
length_element = sqrt(xa*xa+ya*ya);
cosa = xa/length_element;
sena = ya/length_element;
ll = length_element;
L = [cosa*eye(2) sena*eye(2) zeros(2);
-sena*eye(2) cosa*eye(2) zeros(2);
zeros(2,4) eye(2)];
oneu = 1/6*[2 1;1 2];
oneu2 = 1/70*[26 9; 9 26];
oneu3 = ll/420*[22 -13; 13 -22];
oneu4 = ll^2/420*[4 -3; -3 4];
m1 = rho_M*ll*[oneu zeros(2,4);
zeros(2) oneu2 oneu3;
zeros(2) oneu3' oneu4];
mass(elementDof,elementDof) = mass(elementDof,elementDof) + L'*m1*L;
end

% Droppers

for e = 2*numberElements_W+1:(2*numberElements_W+9*No_span)
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [indice indice+Nodes_effective indice+2*Nodes_effective];
    xa = xx(indice(2))-xx(indice(1));
ya = yy(indice(2))-yy(indice(1));
length_element = sqrt(xa*xa+ya*ya);
cosa = xa/length_element;
sena = ya/length_element;
ll = length_element;
L = [cosa*eye(2) sena*eye(2) zeros(2);
-sena*eye(2) cosa*eye(2) zeros(2);
zeros(2,4) eye(2)];
oneu = 1/6*[2 1;1 2];
oneu2 = 1/70*[26 9; 9 26];
oneu3 = ll/420*[22 -13; 13 -22];
oneu4 = ll^2/420*[4 -3; -3 4];
m1 = rho_D*ll*[oneu zeros(2,4);
zeros(2) oneu2 oneu3;
zeros(2) oneu3' oneu4];
mass(elementDof,elementDof) = mass(elementDof,elementDof) + L'*m1*L;
end
% 240, 19
% 465, 19
% 911, 19
% 1132, 19

%% Prescribed Degree of freedom
% Prescribed Degree of freedom for supports
prescribedDof=zeros(1,6*(2*No_span+2)); i=0;
for e =(2*numberElements_W+9*No_span)+1:(No_span*9+2*numberElements_W)+2+2*No_span
       indice = elementNodes(e,:);
        prescribedDof1 = [indice indice+Nodes_effective indice+2*Nodes_effective];
              for j=1:6
               prescribedDof(1,j+i)=prescribedDof1(1,j);
              end
           i=i+6;
end
% Prescribed Degree of freedom for droppers
for e=(2*numberElements_W+1):(2*numberElements_W+No_droppers*No_span)
    indice=elementNodes(e,:);
    prescribedDof2(e-2*numberElements_W)=[indice(1)+Nodes_effective];
end
% nullifying matrices to avoid sigularity error
   for i=0:5:(length(prescribedDof)-(2*No_span+4))
      prescribedDof(:,3+i)=[];
   end
% Prescribed Degree of freedom for droppers  
  q=length(prescribedDof);
   for i=length(prescribedDof)+1:length(prescribedDof)+No_droppers*No_span
       prescribedDof(i)=prescribedDof2(i-q);
   end

%% Solution

% solution

activeDof = setdiff((1:GDof)', prescribedDof);
s_activ = stiffness(activeDof, activeDof);
f_activ = force(activeDof);
U = stiffness(activeDof,activeDof)\force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;

% reactions
F1 = stiffness*displacements;
reactions = F1(prescribedDof);

[prescribedDof' reactions];

q=displacements(((Nodes_effective+numberNodes_W)+1):(((Nodes_effective+numberNodes_W)+1)+numberElements_W));
for i=1:numberNodes_W
    MW_selfweight(i)=NC_MW(i,2)+q(i);
end

%% Before applying self weight of the wires
figure()
ln_c=plot(NC_CW(:,1),NC_CW(:,2),'k','LineWidth',4) % Plot for the contact wire
ln_c.LineWidth = 4;
ln_c.Color = [0.6350 0.0780 0.1840];
hold on
ln_m=plot(NC_MW(:,1),NC_MW(:,2),'m','LineWidth',4) % Plot for the messenger wire
ln_m.LineWidth = 4;
ln_m.Color = [0.4660 0.6740 0.1880];

% plot for the droppers
for i=1:2:(length(NC_droppers)-1)
    plot([NC_droppers(i,1) NC_droppers(i+1,1)],[NC_droppers(i,2) NC_droppers(i+1,2)],'r','LineWidth',2)
end 
xlabel('Length of wire [m]')
ylabel('Vertical Position [m]')
xlim([-20 No_span*lengthSpan+20])
ylim([-0.2 1.4])
hold off

%% Plot for the dropper length and dropper force
figure()
bar(F,0.6)
xlabel('Droppers')
ylabel('Force [N]')
title('Dropper force')

figure()
bar(length_dropper,0.6)
xlabel('Droppers')
ylabel('Length of dropper [m]')
title('Length of the droppers')

%% Plotting the results for the static condition after applying self weight
figure()

% plot for contact wire

% adding the initial configuration of contact wire to the results obtained
for j=1:No_span
    for i=1:221
        if j==1
            displacements(Nodes_effective+i)=displacements(Nodes_effective+i)+CW_Y(i);
            
        else
            displacements((Nodes_effective+i+(j-1)*221)-(j-1))= displacements((Nodes_effective+i+(j-1)*221)-(j-1))+CW_Y(i);
            
        end
    end
end
   
ln_c=plot(xx(1:numberNodes_W),displacements((Nodes_effective+1):((Nodes_effective+1)+numberElements_W)))
ln_c.LineWidth = 2;
ln_c.Color = [0.6350 0.0780 0.1840];

% Plot for the messenger wire
hold on
ln_m=plot(xx(1:numberNodes_W),MW_selfweight)
ln_m.LineWidth = 2;
ln_m.Color = [0.4660 0.6740 0.1880];
hold off
xlabel('Length of wire [m]')
ylabel('Vertical Position [m]')
xlim([-20 No_span*lengthSpan+20])
ylim([-0.2 1.4])

title("Static configuration of contact and messenger wire for " + No_span + " spans ")