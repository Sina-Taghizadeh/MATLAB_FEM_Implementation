% A general MATLAB code to solve 2D problems with truss, triangle and bilinear quadrilateral elements
% by Sina Taghizadeh
% July 2022
author = "Sina Taghizadeh";
license = "GNU GPL Version 3.0";

clc;clear;close all

%----------------------------------------------------------------------
% User Input or Automatic Mesh Generation for Triangle Elements
%----------------------------------------------------------------------

% insert the excel file name
Data='input_sample_truss.xlsx';
excelfile=Data;

% X,y of each node in global cordinate system
Nodes_xy=xlsread(excelfile,'k2:oo3')'; 

% each row is these informations of each element: 
% [ first node number, second node number, third node number, 4th node number]
elements=xlsread(excelfile,'b2:e1000');

% other properties of each element
TrussElementsArea=xlsread(excelfile,'f2:f1000');
EandNu=xlsread(excelfile,'h2:i1000');
ElementsThickness=xlsread(excelfile,'g2:g1000');

% Given all external forces 
Given_Force=xlsread(excelfile,'k4:oo5');

% Given Fixed Boundary Conditions 
zero_dof=xlsread(excelfile,'k6:oo7');  

% Automatic mesh generation just for triangle elements and rectangular domains if needed
generation=0;   % should fill by user (1 for on & 0 for off) 
if generation
    RectanguLarlength=9;   %The length of the rectangular domain
    RectanguLarWidth=1;    %The width of the rectangular domain
    n_length=20;   %should fill by user
    n_width=8;     %should fill by user
    Nodes_gen=zeros((n_length+1)*(n_width/2+1),2);  %preallocation
    elem_gen=zeros(n_length*n_width,3);             %preallocation
    a=n_width/2+1;     %Auxiliary variable
    b=0;             %Counter
    for i=1:a:size(Nodes_gen)      %for X coordinate of nodes
        Nodes_gen(i,1)=b*RectanguLarlength/n_length;
        Nodes_gen(i+1:i+a-1,1)=Nodes_gen(i,1);
        b=b+1;
    end
    b=0;
    for i=1:a                      %for y coordinate of nodes
        Nodes_gen(i,2)=2*b*RectanguLarWidth/n_width;
        b=b+1;
        for o=1:n_length
            Nodes_gen(i+o*a,2)=Nodes_gen(i,2);
        end
    end
    b=1;
    for i=1:2:n_width*n_length         %for node connectivity of elements
        elem_gen(i:i+1,1)=b;
        elem_gen(i+1,3)=elem_gen(i+1,1)+1;
        elem_gen(i+1,2)=b+a+1;
        elem_gen(i,2)=elem_gen(i+1,1)+a;
        elem_gen(i,3)=b+a+1;
        if mod(b,a)==n_width/2
            b= b+1;
        end
        b=b+1;
    end
    Nodes_xy=Nodes_gen;
    elements=elem_gen;
    
    % Change the rest of the properties of the Excel file
    EandNu(3:n_length*n_width,1)=EandNu(1,1);
    EandNu(3:n_length*n_width,2)=EandNu(1,2);
    ElementsThickness(3:n_length*n_width)=ElementsThickness(1);
    Given_Force(:,1:(n_length+1)*(n_width/2+1))=0;
    Given_Force(1,(n_length+1)*(n_width/2+1))=+100000;
    Given_Force(1,(n_length+1)*(n_width/2+1)-n_width/2)=-100000;
    zero_dof(:,1:a)=1;
    zero_dof(:,a+1:(n_length+1)*(n_width/2+1))=0;
end
%--------------------------------------------------------------------------
% End of User Input or Automatic Mesh Generation for Triangle Elements
%--------------------------------------------------------------------------


%----------------------------------------------------------------------
% Production of Problem Prerequisites
%----------------------------------------------------------------------

tic;   %for runtime of my program

null = 0.000001; %define zero for problem
ScalingFactor=100; %for Visualization of results

% Number of nodes in the problem
n_nodes=size(Nodes_xy,1);

% total degree of fredom (for planar problems)
tot_dof=n_nodes * 2;

% number of elements
no_elements = size(elements,1);

% preallocation for increase speed of program in loops
displacements= zeros(tot_dof,1);   % displacement vector
stiffness=zeros(tot_dof);          % global stiffness matrix
Force=zeros(no_elements,1);       % stress vector

% P = matrix of each node cordinate numbers in Global
P=reshape(1:tot_dof,2,[])';
%--------------------------------------------------------------------------
% End of Production of Problem Prerequisites
%--------------------------------------------------------------------------


%----------------------------------------------------------------------
% Stiffness Matrix Generation and Assembly
%----------------------------------------------------------------------

% creat stiffnes matrices and assembly for 2D elements
for     j=1:no_elements    %loop for creating k element
    
    if size(elements,2) == 2  %truss elements
        L= ((Nodes_xy(elements(j,2),1) - Nodes_xy(elements(j,1),1))^2 + (Nodes_xy(elements(j,2),2) - Nodes_xy(elements(j,1),2))^2  )^0.5 ;
        c=(Nodes_xy(elements(j,2),1) - Nodes_xy(elements(j,1),1))/L;         %cos    
        s=(Nodes_xy(elements(j,2),2) - Nodes_xy(elements(j,1),2))/L;         %sin     
        Kelement=(TrussElementsArea(j)*EandNu(j)/L)*[c*c c*s -c*c -c*s;c*s s*s -c*s -s*s;-c*c -c*s c*c c*s;-c*s -s*s c*s s*s]; % k for each element
        stiffness([P(elements(j,1),:) P(elements(j,2),:)],[P(elements(j,1),:) P(elements(j,2),:)])= stiffness([P(elements(j,1),:) P(elements(j,2),:)],[P(elements(j,1),:) P(elements(j,2),:)]) +  Kelement; %assembly
    end
    
    if size(elements,2) == 3  %triangle elements
        delta=0.5*det([1 Nodes_xy(elements(j,1),1) Nodes_xy(elements(j,1),2); 1 Nodes_xy(elements(j,2),1) Nodes_xy(elements(j,2),2); 1 Nodes_xy(elements(j,3),1) Nodes_xy(elements(j,3),2)]);
        Belement=(1/(2*delta))*[Nodes_xy(elements(j,2),2)-Nodes_xy(elements(j,3),2) 0 Nodes_xy(elements(j,3),2)-Nodes_xy(elements(j,1),2) 0 Nodes_xy(elements(j,1),2)-Nodes_xy(elements(j,2),2) 0 ;
                                0 Nodes_xy(elements(j,3),1)-Nodes_xy(elements(j,2),1) 0 Nodes_xy(elements(j,1),1)-Nodes_xy(elements(j,3),1) 0 Nodes_xy(elements(j,2),1)-Nodes_xy(elements(j,1),1) ;
                                Nodes_xy(elements(j,3),1)-Nodes_xy(elements(j,2),1) Nodes_xy(elements(j,2),2)-Nodes_xy(elements(j,3),2) Nodes_xy(elements(j,1),1)-Nodes_xy(elements(j,3),1) Nodes_xy(elements(j,3),2)-Nodes_xy(elements(j,1),2) Nodes_xy(elements(j,2),1)-Nodes_xy(elements(j,1),1) Nodes_xy(elements(j,1),2)-Nodes_xy(elements(j,2),2)];
        C=(EandNu(j,1)/(1-(EandNu(j,2))^2))*[1 EandNu(j,2) 0; EandNu(j,2) 1 0; 0 0 (1-EandNu(j,2))/2]; %for plane stress
        %C=(EandNu(j,1)/((1+(EandNu(j,2)))*(1-2*(EandNu(j,2)))))*[1-EandNu(j,2) EandNu(j,2) 0; EandNu(j,2) 1-EandNu(j,2) 0; 0 0 (1-2*(EandNu(j,2)))/2]; %for plane strain
        Kelement=delta*ElementsThickness(j)*transpose(Belement)*C*Belement;
        stiffness([P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:)],[P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:)])= stiffness([P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:)],[P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:)]) +  Kelement; %assembly
    end 
    
    if size(elements,2) == 4  %bilinear quadrilateral elements
        syms e n
        N=[1/4*(1-e)*(1-n);
            1/4*(1+e)*(1-n);
            1/4*(1+e)*(1+n);
            1/4*(1-e)*(1+n)];
        x=N(1)*Nodes_xy(elements(j,1),1)+N(2)*Nodes_xy(elements(j,2),1)+N(3)*Nodes_xy(elements(j,3),1)+N(4)*Nodes_xy(elements(j,4),1);
        y=N(1)*Nodes_xy(elements(j,1),2)+N(2)*Nodes_xy(elements(j,2),2)+N(3)*Nodes_xy(elements(j,3),2)+N(4)*Nodes_xy(elements(j,4),2);
        jac=[diff(x,e) diff(y,e);diff(x,n) diff(y,n)];
        for o=1:4
            a=inv(jac)*[diff(N(o),e);diff(N(o),n)];
            bder(o)=a(1);
            bder(o+4)=a(2);
        end
        Belement=[bder(1) 0 bder(2) 0 bder(3) 0 bder(4) 0; 0 bder(5) 0 bder(6) 0 bder(7) 0 bder(8); bder(5) bder(1) bder(6) bder(2) bder(7) bder(3) bder(8) bder(4) ];
        C=(EandNu(j,1)/(1-(EandNu(j,2))^2))*[1 EandNu(j,2) 0; EandNu(j,2) 1 0; 0 0 (1-EandNu(j,2))/2]; %for plane stress
        %C=(EandNu(j,1)/((1+(EandNu(j,2)))*(1-2*(EandNu(j,2)))))*[1-EandNu(j,2) EandNu(j,2) 0; EandNu(j,2) 1-EandNu(j,2) 0; 0 0 (1-2*(EandNu(j,2)))/2]; %for plane strain
        Kelement=0;
        for t=1:4
            gauss=[0.57735 0.57735;-0.57735 0.57735;0.57735 -0.57735;-0.57735 -0.57735];
            n=gauss(t,1);
            e=gauss(t,2);
            knew=subs(transpose(Belement)*C*Belement*det(jac));
            Kelement = Kelement+ knew;
        end
        stiffness([P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:) P(elements(j,4),:)],[P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:) P(elements(j,4),:)])= stiffness([P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:) P(elements(j,4),:)],[P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:) P(elements(j,4),:)]) +  Kelement; %assembly
    end
end
%--------------------------------------------------------------------------
% End of Stiffness Matrix Generation and Assembly
%--------------------------------------------------------------------------


%----------------------------------------------------------------------
% Apply BCs
%----------------------------------------------------------------------

% apply zero dof to given dof for solve
Given_dof=find(zero_dof(:)');

% nodes with unkhown field variables
unknown_dofNodes= setdiff((1:tot_dof)',Given_dof');

% apply given external forces to force vector
Forces = Given_Force(:)  ;
%--------------------------------------------------------------------------
% End of Apply BCs
%--------------------------------------------------------------------------


%----------------------------------------------------------------------
% Solve Linear System with LU Factorization Decomposition
%----------------------------------------------------------------------

% determination of unkhown displacements with LU decomposition
[Li,Ui,qi] = lu(stiffness(unknown_dofNodes,unknown_dofNodes));
y = Li\(qi*Forces(unknown_dofNodes));
displacements(unknown_dofNodes) = Ui\y;
u=displacements;
%--------------------------------------------------------------------------
% End of Solve Linear System with LU Factorization Decomposition
%--------------------------------------------------------------------------


%----------------------------------------------------------------------
% Post-Processing and Visualization of Results
%----------------------------------------------------------------------

% determination of nodal forces
F= stiffness * displacements;
F(abs(F)<null)=0;

% determination of Forces and stresse of each elements for 2D
if size(elements,2) == 2  %truss element
    for     j=1:no_elements
    L= ((Nodes_xy(elements(j,2),1) - Nodes_xy(elements(j,1),1))^2 + (Nodes_xy(elements(j,2),2) - Nodes_xy(elements(j,1),2))^2  )^0.5 ;
    c=(Nodes_xy(elements(j,2),1) - Nodes_xy(elements(j,1),1))/L; %cx
    s=(Nodes_xy(elements(j,2),2) - Nodes_xy(elements(j,1),2))/L; %cy
    Force(j)=(EandNu(j)/L)*[-c -s c s]*u([P(elements(j,1),:) P(elements(j,2),:)])*TrussElementsArea(j);
    Force(abs(Force)<null)=0;
    end
end
if size(elements,2) == 3  %triangle elements
    sigma=zeros(no_elements,3); %preallocation
    for     j=1:no_elements
    delta=0.5*det([1 Nodes_xy(elements(j,1),1) Nodes_xy(elements(j,1),2); 1 Nodes_xy(elements(j,2),1) Nodes_xy(elements(j,2),2); 1 Nodes_xy(elements(j,3),1) Nodes_xy(elements(j,3),2)]);
    Belement=(1/(2*delta))*[Nodes_xy(elements(j,2),2)-Nodes_xy(elements(j,3),2) 0 Nodes_xy(elements(j,3),2)-Nodes_xy(elements(j,1),2) 0 Nodes_xy(elements(j,1),2)-Nodes_xy(elements(j,2),2) 0 ;
                                0 Nodes_xy(elements(j,3),1)-Nodes_xy(elements(j,2),1) 0 Nodes_xy(elements(j,1),1)-Nodes_xy(elements(j,3),1) 0 Nodes_xy(elements(j,2),1)-Nodes_xy(elements(j,1),1) ;
                                Nodes_xy(elements(j,3),1)-Nodes_xy(elements(j,2),1) Nodes_xy(elements(j,2),2)-Nodes_xy(elements(j,3),2) Nodes_xy(elements(j,1),1)-Nodes_xy(elements(j,3),1) Nodes_xy(elements(j,3),2)-Nodes_xy(elements(j,1),2) Nodes_xy(elements(j,2),1)-Nodes_xy(elements(j,1),1) Nodes_xy(elements(j,1),2)-Nodes_xy(elements(j,2),2)];
    C=(EandNu(j,1)/(1-(EandNu(j,2))^2))*[1 EandNu(j,2) 0; EandNu(j,2) 1 0; 0 0 (1-EandNu(j,2))/2];
    %C=(EandNu(j,1)/((1+(EandNu(j,2)))*(1-2*(EandNu(j,2)))))*[1-EandNu(j,2) EandNu(j,2) 0; EandNu(j,2) 1-EandNu(j,2) 0; 0 0 (1-2*(EandNu(j,2)))/2];
    sigma(j,:)=C*Belement*u([P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:)]);
    sigma(abs(sigma)<null)=0;
    end
end
if size(elements,2) == 4  %bilinear quadrilateral elements
    sigma=zeros(no_elements,3); %preallocation
    for     j=1:no_elements
        syms e n
        N=[1/4*(1-e)*(1-n);
            1/4*(1+e)*(1-n);
            1/4*(1+e)*(1+n);
            1/4*(1-e)*(1+n)];
        x=N(1)*Nodes_xy(elements(j,1),1)+N(2)*Nodes_xy(elements(j,2),1)+N(3)*Nodes_xy(elements(j,3),1)+N(4)*Nodes_xy(elements(j,4),1);
        y=N(1)*Nodes_xy(elements(j,1),2)+N(2)*Nodes_xy(elements(j,2),2)+N(3)*Nodes_xy(elements(j,3),2)+N(4)*Nodes_xy(elements(j,4),2);
        jac=[diff(x,e) diff(y,e);diff(x,n) diff(y,n)];
        for o=1:4
            a=inv(jac)*[diff(N(o),e);diff(N(o),n)];
            bder(o)=a(1);
            bder(o+4)=a(2);
        end
        Belement=[bder(1) 0 bder(2) 0 bder(3) 0 bder(4) 0; 0 bder(5) 0 bder(6) 0 bder(7) 0 bder(8); bder(5) bder(1) bder(6) bder(2) bder(7) bder(3) bder(8) bder(4) ];
        C=(EandNu(j,1)/(1-(EandNu(j,2))^2))*[1 EandNu(j,2) 0; EandNu(j,2) 1 0; 0 0 (1-EandNu(j,2))/2];
        w=zeros(4,3);
        for t=1:4
            gauss=[0.57735 0.57735;-0.57735 0.57735;0.57735 -0.57735;-0.57735 -0.57735];
            n=gauss(t,1);
            e=gauss(t,2);
            CB=subs(C*Belement);
            w(t,:)=CB*u([P(elements(j,1),:) P(elements(j,2),:) P(elements(j,3),:) P(elements(j,4),:)]); %stress in each gauss point
        end
        sigma(j,:)=max(w); %find max stress between stresses in gauss points
        sigma(abs(sigma)<null)=0;
    end
end
toc;

% Visualization of results for truss elements
if size(elements,2) == 2  %truss elements
    disp('all of displacemnts and reaction forces of this problem:');
    
    %find coordinate of nodes after deformation for plot
    v=Nodes_xy+ScalingFactor*reshape(u,2,[])';  %v is node coordinates after deformation
    
    nodes=zeros(2, n_nodes);
    for a=1:n_nodes
        nodes(:,a)=a;
    end
    nodes=nodes(:);
    orientation=repmat(['u';'v'],n_nodes,1);
 
    
    T=table(orientation,nodes,u,F)
    disp('element stresses');
    Elements=(1:no_elements)';
    t=table(Elements,Force)
    
    % plot trust before and after deformation 

    %plot trust before deformation
    for w=1:no_elements
    p1=plot(Nodes_xy(elements(w,1:2),1),Nodes_xy(elements(w,1:2),2),'o-','Color','k','linewidth',2,'markersize',4);
    hold on
    end
    
    hold on
    
    %plot trust after deformation
    for s=1:no_elements
        p2=plot(v(elements(s,1:2),1),v(elements(s,1:2),2),'o--','Color','red','linewidth',2,'markersize',4);
        hold on
    end
    
    legend([p1 p2],{'before',"after (Scaling Factor = " + ScalingFactor + ")"})
 
end

% Visualization of results for triangle elements
if size(elements,2) == 3  %triangle elements
    disp('all of displacemnts and stresses of this problem:(for CST element)');
    
    %find coordinate of nodes after deformation for plot
    v=Nodes_xy+ScalingFactor*reshape(u,2,[])';  %v is node coordinates after deformation
    
    nodes=zeros(2, n_nodes);
    for a=1:n_nodes
        nodes(:,a)=a;
    end
    nodes=nodes(:);
    orientation=repmat(['u';'v'],n_nodes,1);
 
    
    T=table(orientation,nodes,u,F)
    disp('element stresses');
    Elements=(1:no_elements)';
    t=table(Elements,sigma)
    
    % plot before and after deformation 

    %plot before deformation
    for w=1:no_elements
        p1=plot(Nodes_xy(elements(w,1:2),1),Nodes_xy(elements(w,1:2),2),'o-','Color','k','linewidth',2,'markersize',4);
        hold on
        p2=plot(Nodes_xy(elements(w,2:3),1),Nodes_xy(elements(w,2:3),2),'o-','Color','k','linewidth',2,'markersize',4);
        hold on
        p3=plot(Nodes_xy([elements(w,1) elements(w,3)],1),Nodes_xy([elements(w,1) elements(w,3)],2),'o-','Color','k','linewidth',2,'markersize',4);
    end
    
    hold on
    
    %plot after deformation
    for s=1:no_elements
        p4=plot(v(elements(s,1:2),1),v(elements(s,1:2),2),'o--','Color','red','linewidth',2,'markersize',4);
        hold on
        p5=plot(v(elements(s,2:3),1),v(elements(s,2:3),2),'o--','Color','red','linewidth',2,'markersize',4);
        hold on
        p6=plot(v([elements(s,1) elements(s,3)],1),v([elements(s,1) elements(s,3)],2),'o--','Color','red','linewidth',2,'markersize',4);
    end
    
    legend([p1 p4],{'before',"after (Scaling Factor = " + ScalingFactor + ")"})
 
end

% Visualization of results for bilinear quadrilateral elements
if size(elements,2) == 4   %bilinear quadrilateral elements
    disp('all of displacemnts and stresses of this problem:(for BILINEAR QUADRILATERAL element)');
    
    %find coordinate of nodes after deformation for plot
    v=Nodes_xy+ScalingFactor*reshape(u,2,[])';  %v is node coordinates after deformation
    
    nodes=zeros(2, n_nodes);
    for a=1:n_nodes
        nodes(:,a)=a;
    end
    nodes=nodes(:);
    orientation=repmat(['u';'v'],n_nodes,1);
 
    
    T=table(orientation,nodes,u,F)
    disp('element stresses');
    Elements=(1:no_elements)';
    t=table(Elements,sigma)
    
    % plot before and after deformation 

    %plot before deformation
    for w=1:no_elements
        p1=plot(Nodes_xy(elements(w,1:2),1),Nodes_xy(elements(w,1:2),2),'o-','Color','k','linewidth',2,'markersize',4);
        hold on
        p2=plot(Nodes_xy(elements(w,2:3),1),Nodes_xy(elements(w,2:3),2),'o-','Color','k','linewidth',2,'markersize',4);
        hold on
        p3=plot(Nodes_xy(elements(w,3:4),1),Nodes_xy(elements(w,3:4),2),'o-','Color','k','linewidth',2,'markersize',4);
        hold on
        p4=plot(Nodes_xy([elements(w,1) elements(w,4)],1),Nodes_xy([elements(w,1) elements(w,4)],2),'o-','Color','k','linewidth',2,'markersize',4);
    end
    
    hold on
    
    %plot after deformation
    for s=1:no_elements
        p5=plot(v(elements(s,1:2),1),v(elements(s,1:2),2),'o--','Color','red','linewidth',2,'markersize',4);
        hold on
        p6=plot(v(elements(s,2:3),1),v(elements(s,2:3),2),'o--','Color','red','linewidth',2,'markersize',4);
        hold on
        p7=plot(v(elements(s,3:4),1),v(elements(s,3:4),2),'o--','Color','red','linewidth',2,'markersize',4);
        hold on
        p8=plot(v([elements(s,1) elements(s,4)],1),v([elements(s,1) elements(s,4)],2),'o--','Color','red','linewidth',2,'markersize',4);
    end
    
    legend([p1 p5],{'before',"after (Scaling Factor = " + ScalingFactor + ")"})
 
end
%--------------------------------------------------------------------------
% End of Post-Processing and Visualization of Results
%--------------------------------------------------------------------------
