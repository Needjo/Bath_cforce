% Configurational forces 2d problem - Bath 27.10.2014
% Crack propagation by remeshing using configurational force analysis. 
% Example problem from Miehe C., Gurses E., "A robust algorithm for 
% configurational-force-driven brittle crack propagation with R-adaptive 
% mesh alignment"-example 4.1.3.Single edge notched specimen under tension.
% This procedure simulates a crack propagation using configurational force
% driven remeshing. The model is purely kinematical.
% Procedure author: Nikola Lustig, mag.ing.aedif.

clear; 
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters

lgrede = 2;
hgrede = 0.5;
tgrede = 1;

nlelemenata = 40;
nhelemenata = 10;

notch = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notch_direction describes the notch position direction
% 1-horizontal notch propagating from the left edge
% 2-horizontal notch propagating from the right edge
% 3-vertical notch propagating from the lower edge
% 4-vertical notch propagating from the upper edge
Notch_direction = 1; 


E = 2 * 80000000000 * ( 1 + 0.3 );
Nu = 0.3;

% Lame constants calculation for comparison with paper
G = E/(2*(1+Nu));
Lame_1 = E * Nu / ( ( 1 + Nu ) * ( 1 - 2 * Nu ) );

% p0 defines a distributed load
p0 = 10^7;
% Sila defines a concentrated load
% Sila = -100;

max_koraka = 16;
tol = 10^-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite element mesh import

textfilename = ['Mesh files\Q4_mesh_SEN_extended_' num2str(lgrede) '_' num2str(hgrede) '_' num2str(notch) '_EL_' num2str(nlelemenata) 'x' num2str(nhelemenata) '_DIR_' num2str(Notch_direction) '.txt'];
[ xy , EL , nC , nEL , nSEG , SEG , EL_SEG , EL_NEIGHBOUR , SEG_NEIGHBOUR , Node_occurence_EL , Node_occurence_SEG , Node_EL , Node_SEG ] = Q4_mesh_reader_extended ( textfilename );

ELX = zeros ( 5 , nEL );
ELY = zeros ( 5 , nEL );
for i = 1:nEL
    ELX ( 1 , i ) = EL ( 1 , i );
    ELX ( 2 , i ) = xy ( 1, EL ( 2 , i ) );
    ELX ( 3 , i ) = xy ( 1, EL ( 3 , i ) );
    ELX ( 5 , i ) = xy ( 1, EL ( 5 , i ) );
    ELX ( 4 , i ) = xy ( 1, EL ( 4 , i ) );
end;
for i = 1:nEL
    ELY ( 1 , i ) = EL ( 1 , i );
    ELY ( 2 , i ) = xy ( 2, EL ( 2 , i ) );
    ELY ( 3 , i ) = xy ( 2, EL ( 3 , i ) );
    ELY ( 5 , i ) = xy ( 2, EL ( 5 , i ) );
    ELY ( 4 , i ) = xy ( 2, EL ( 4 , i ) );
end;


ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT = [ELX(2:5,:);ELX(2,:)];
ELYPLOT = [ELY(2:5,:);ELY(2,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plane strain state for isotropic linear elastic material
Smatrica = [ (1-Nu^2)/E -Nu*(1+Nu)/E 0; -Nu*(1+Nu)/E (1-Nu^2)/E 0; 0 0 1/G ];
Cmatrica = inv(Smatrica);

% % Plane stress state for isotropic linear elastic material
% Cmatrica = ( E / ( 1 - Nu^2 ) ) * [ 1 Nu 0; Nu 1 0; 0 0 (1-Nu)/2 ]; 


% Gaussian numerical quadrature rule - 2x2 integration
[ rgauss, sgauss, wgauss ] = gaussquad2_shift;
rlength = length ( rgauss );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary values and loading definition

% Rubni = [ 1 , find(abs(xy(1,:)-lgrede) <= tol & xy(2,:)==0) , find(xy(1,:)==lgrede/2 & xy(2,:)==hgrede) ];
Rubni = [ find(abs(xy(1,:) <= tol) & xy(2,:)==hgrede) , find(abs(xy(1,:)-lgrede) <= tol & xy(2,:)==hgrede) , find(xy(1,:)==lgrede/2 & xy(2,:)==0) ];
Rubni = [ Rubni ; 1 1 0 ; 0 0 1 ];

% Use this for DCB
Rubni = [ find(abs(xy(1,:) - lgrede ) <= tol) ];
Rubni = [ Rubni ; zeros(2,length(Rubni)) ];


FV = zeros ( 2 * nC , 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a vertical end-node concentrated load
% FV ( 2 * nC ) = Sila;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a DCB
Sila = p0;
% FV ( 2 * nC ) = Sila;
% FV ( 2 * ( nlelemenata + 1 ) ) = -Sila;
FV ( 2 ) = -Sila;
FV ( 2*(nC - nlelemenata) ) = Sila;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a vertical distributed load
% Sila = p0 / (nlelemenata+1);
% FV ( 2*(nC - (1:nlelemenata+1)) ) = Sila;
% FV(2*nC) = Sila/2;
% FV(2*(nC-nlelemenata))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a horizontal distributed load on both sides
% Sila = p0 * hgrede / (nhelemenata);
% 
% FV ( 2*find(xy(1,:)==round(lgrede*100)/100) - 1 ) = Sila;
% FV ( 2*find(xy(1,:)==0) - 1 ) = -Sila;
% 
% FV(2*nC - 1) = Sila/2;
% FV(2*(nlelemenata+2)-1)=Sila/2;
% FV(1) = -Sila/2;
% FV(2*(nC-nlelemenata)-1) = -Sila/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation of data for analysis

Naprezanja = zeros ( nEL * rlength , 5 );
Deformacije = zeros ( nEL * rlength , 5 );
Def_grad_lista = zeros ( nEL * rlength , 4 );
W = zeros ( nEL * rlength , 1 );
Eshelby2 = zeros ( nEL * rlength , 4 );
shape4lin_save = zeros ( nEL * rlength , 4 );
dxshape4lin_save = zeros ( nEL * rlength , 4 );
dyshape4lin_save = zeros ( nEL * rlength , 4 );
detJacob_save = zeros ( nEL * rlength , 1 );
Bmatrica_save = zeros ( nEL * rlength , 24 );
Bmatrica2_save = zeros ( nEL * rlength , 32 );
Pomaci_elementa = zeros ( 8 , 1 );
xy_popis = cell ( max_koraka , 1 );
pomak_popis = cell ( max_koraka , 1 );
Conf_popis = cell ( max_koraka , 1 );
Energy_FEM_popis = zeros ( 1 , max_koraka );
norma_popis = zeros ( 1 , max_koraka );
J_integral_num = zeros ( 1 , max_koraka );
Oduzmi = 0;

Usporedba_korak = [  ];
Naprezanja_korak = [  ];
Deformacije_korak = [  ];
detJacob_korak = [  ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of loop over elements for mesh improvement - NODE SHIFTING

for broj_koraka = 1 : max_koraka
    
    if broj_koraka ~= 1
        
        [ J_integral , Crack_tip ] = max ( abs ( Conf_force ) );
        
% %         if broj_koraka==2, 
% %             Crack_tip = 2*find(round(xy ( 1 , : )*100)/100 == lgrede/2 & xy ( 2 , : ) == 0);
% % %             Crack_tip = 2*find(round(xy ( 1 , : )*100)/100 == lgrede/2 & round(xy ( 2 , : )*100)/100 == hgrede);
% %         end;
        
        if rem ( Crack_tip , 2 ) == 0
            
            Crack_tip = Crack_tip / 2;
            Crack_propagation_direction = 0;
            
        else
            
            Crack_tip = ( Crack_tip + 1 ) / 2;
            Crack_propagation_direction = 1;
            
        end;
        
        J_integral = -J_integral;
                
%         if broj_koraka == 2
%             Seg_tip = 23; 
%         elseif broj_koraka == 3
%             Seg_tip = 28;
%         elseif broj_koraka == 4
%             Seg_tip = 33;
%         else
%             Seg_tip = 38;
%         end;

        Crack_tip_vector = [ Conf_force(2*Crack_tip-1,1); Conf_force(2*Crack_tip,1) ] .* ( -1 ); % This way we turn the vector by 180 degrees
        Crack_tip_angle = atan ( Conf_force(2*Crack_tip,1) / Conf_force(2*Crack_tip-1,1) );
        J_magnitude = norm( [ Conf_force(2*Crack_tip,1) Conf_force(2*Crack_tip-1,1) ] );
        
        
        Crack_angle_nodes = sum ( SEG ( : , Node_SEG ( Crack_tip , Node_SEG ( Crack_tip , : ) > 0 ) ) ) - Crack_tip;
        
        for ii = 1 : length ( Crack_angle_nodes ) 
            
            Crack_propagation_vector ( : , ii ) = xy ( : , Crack_angle_nodes(ii) ) - xy ( : , Crack_tip );
            Crack_propagation_magnitude ( ii ) = norm ( Crack_propagation_vector ( : , ii ) );
            Crack_propagation_angle ( ii ) = atan ( Crack_propagation_vector ( 2 , ii ) / Crack_propagation_vector ( 1 , ii ) );
            
        end;
       
       for ii = 1 : length ( Crack_angle_nodes )  
            Crack_angles ( ii ) = acos ( sum ( Crack_propagation_vector(:,ii) .* Crack_tip_vector ) ./ ( Crack_propagation_magnitude(ii) .* J_magnitude ) );
       end;
       [ A , Crack_angle_critical ] = min ( abs ( Crack_angles ) ); 
       
%         [ A , Crack_angle_critical ] = min ( abs ( Crack_propagation_angle - Crack_tip_angle - pi ) );
        
        Seg_tip = find(sum(SEG == repmat([ Crack_tip; Crack_angle_nodes(Crack_angle_critical) ],1,nSEG))==2 | sum(SEG == repmat([ Crack_angle_nodes(Crack_angle_critical);Crack_tip ],1,nSEG))==2);
        
% %         if broj_koraka==2, 
% %             Seg_tip = find(sum(repmat([Crack_tip;Crack_tip+nlelemenata+1],1,nSEG)==SEG)==2);
% % %             Seg_tip = find(sum(repmat([Crack_tip-nlelemenata-1;Crack_tip],1,nSEG)==SEG)==2);
% %         end;
 
%         A SEGMENT SELECTION TECHNIQUE IS NEEDED!
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Starting the crack propagation by invoking segment and node
        % doubling and appropriate data structure updates. Crack_tip node
        % is the node to be doubled while Seg_tip is the segment along
        % which the crack propagation occurs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining the elements connected by the crack tip
        Element_tip_1 = SEG_NEIGHBOUR(1,Seg_tip);
        Element_tip_2 = SEG_NEIGHBOUR(2,Seg_tip);
        
        % Defining the new segment, positioning it at the last place in the 
        % segment list), increasing the global number of segments
        Seg_new = nSEG + 1;
        SEG = [ SEG , SEG(:,Seg_tip) ];
        nSEG = nSEG + 1;
        
        % Changing the neighbouring elements of the new and critical segment
        SEG_NEIGHBOUR ( : , Seg_new ) = [ SEG_NEIGHBOUR(2,Seg_tip) 0 ];
        SEG_NEIGHBOUR ( 2 , Seg_tip ) = 0;

        % Changing the element segment connectivity - the right element is now
        % connected to the new segment. The left element remains as it was.
        % Changing the neighbouring elements for the two elements at the crack tip
        for ii = 1 : 4
    
            if EL_SEG ( ii , Element_tip_2 ) == Seg_tip
        
                EL_SEG ( ii , Element_tip_2 ) = Seg_new;
        
            end;
            
            if EL_NEIGHBOUR ( ii , Element_tip_1 ) == Element_tip_2
        
                EL_NEIGHBOUR ( ii , Element_tip_1 ) = 0;
        
            end;
    
            if EL_NEIGHBOUR ( ii , Element_tip_2 ) == Element_tip_1
        
                EL_NEIGHBOUR ( ii , Element_tip_2 ) = 0;
        
            end;
    
        end;
        
        % Node doubling, the new node is added next to the old node to
        % retain a better global stiffness matrix structure. Increasing the
        % total number of nodes.
        xy = [ xy(:,1:Crack_tip) , xy(:,Crack_tip:nC) ];
        nC = nC + 1;
        New_node = Crack_tip+1;
        
        
        % Renumbering of the nodes due to the addition of a new node in the
        % mesh. Updating the boundary conditions and force vector
        Cracker = EL>Crack_tip;
        Cracker ( 1 , : ) = false;
        EL(Cracker)= EL(Cracker) + 1;
        SEG(SEG>Crack_tip)= SEG(SEG>Crack_tip) + 1;
        for ii = 1:length(Rubni)
            if Rubni(1,ii) > Crack_tip
                Rubni(1,ii) = Rubni(1,ii) + 1;
            end;    
        end;
        FV = [ FV(1:2*Crack_tip,1); [0;0]; FV(2*Crack_tip+1:end) ]; 

        % Finding the element which retains the old node Crack_tip
        Keeper = find ( sum( repmat(Seg_tip,4,nEL) ==  EL_SEG ) > 0 & sum ( repmat(Crack_tip,4,nEL) ==  EL(2:5,:) ) );

        Keeper = [ Keeper EL(1,sum(EL_NEIGHBOUR == Keeper) & sum(EL(2:5,:) == Crack_tip))];
                

        % Logical operators for changing the element topology of the 
        % elements which contain the new node New_node
        logical_aid = true ( 4 , nEL );
        logical_aid ( : , Keeper ) = false;
        logical_aid = repmat(Crack_tip,4,nEL) ==  EL(2:5,:) & logical_aid;
        logical_aid = [ false(1,nEL) ; logical_aid ];
        
        EL(logical_aid) = New_node;

%         Changer = Node_EL ( Crack_tip , Node_EL ( Crack_tip , : ) ~= Keeper & Node_EL ( Crack_tip , : ) ~= 0 );

%         EL ( 2:5 , find ( EL( 2:5 , : ) == Changer ) = New_node;
        
        % Finding the segments which get the new node New_node        
        NonKeeper = find ( sum( repmat(Seg_new,4,nEL) ==  EL_SEG ) > 0 & sum ( repmat(New_node,4,nEL) ==  EL(2:5,:) ) );
%         Node_SEG(Crack_tip,:)~= NonKeeper 
%         NonKeeper = [ NonKeeper EL_SEG(:,SEG_NEIGHBOUR(:,NonKeeper)~=)     EL(1,sum(EL_SEG == NonKeeper) & sum(EL(2:5,:) == Crack_tip))];



        % Logical operators for changing the segment topology of the 
        % segments which contain the new node New_node
        logical_aid = SEG_NEIGHBOUR == NonKeeper;
        logical_aid = (SEG==Crack_tip) & (repmat(logical_aid(1,:) | logical_aid(2,:) ,2,1));

        SEG(logical_aid) = New_node;
        
        
        for ii = 1 : 2
            
            if SEG_NEIGHBOUR ( ii , NonKeeper ) ~= Element_tip_1 & SEG_NEIGHBOUR ( ii , NonKeeper ) ~= Element_tip_2
                
                NonKeeper_element = SEG_NEIGHBOUR ( ii , NonKeeper );
                
            end;
            
        end;
        
        for ii = 1 : 4
            
            if NonKeeper_element == 0
                
                break
                
            end;
            
            for jj = 1 : 2
            
                if SEG(jj,EL_SEG ( ii , NonKeeper_element ) ) == Crack_tip & EL_SEG ( ii , NonKeeper_element ) ~= NonKeeper
                    
                    SEG ( jj,EL_SEG ( ii , NonKeeper_element ) ) = New_node;
                    
                end;
            
            end;
            
        end;

        % Defining the node occurence in elements and segments - regular quadrilateral
        % elements

        Node_occurence_EL = zeros ( nC , 1 );
        Node_occurence_SEG = zeros ( nC , 1 );

        for ii = 1 : nC
            
            Node_occurence_EL ( ii , 1 ) = sum ( sum ( ii == EL ( 2:5 , : ) ) );
            Node_occurence_SEG ( ii , 1 ) = sum ( sum ( ii == SEG  ) );

        end;
        
        % Defining the elements and segments which contain the node
        Node_EL = zeros ( nC , max ( Node_occurence_EL ) );
        Node_SEG = zeros ( nC , max ( Node_occurence_SEG ) );

        for ii = 1 : nC
            
            temp = find ( sum ( ii == EL(2:5,:) ) > 0  );
            
            for jj = 1 : length ( temp )
                
                Node_EL ( ii , jj ) = temp ( jj );
                
            end;
            
            temp = find ( sum ( ii == SEG ) > 0  );
            
            for jj = 1 : length ( temp )
                
                Node_SEG ( ii , jj ) = temp ( jj );
                
            end;
            
        end;

    end;
    
% The new nodes need to be fed into the stiffness matrix    
for ii = 1:nEL
    ELX ( 1 , ii ) = EL ( 1 , ii );
    ELX ( 2 , ii ) = xy ( 1, EL ( 2 , ii ) );
    ELX ( 3 , ii ) = xy ( 1, EL ( 3 , ii ) );
    ELX ( 5 , ii ) = xy ( 1, EL ( 5 , ii ) );
    ELX ( 4 , ii ) = xy ( 1, EL ( 4 , ii ) );
end;
for ii = 1:nEL
    ELY ( 1 , ii ) = EL ( 1 , ii );
    ELY ( 2 , ii ) = xy ( 2, EL ( 2 , ii ) );
    ELY ( 3 , ii ) = xy ( 2, EL ( 3 , ii ) );
    ELY ( 5 , ii ) = xy ( 2, EL ( 5 , ii ) );
    ELY ( 4 , ii ) = xy ( 2, EL ( 4 , ii ) );
end;


ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local stiffness matrix calculation - linear elastic isotropic elements
% Shape functions, their derivatives and Bmatrices saved for later use

Kl = zeros ( 8 , 8 , nEL );

for jj=1:nEL
    
    Klok = zeros ( 8 , 8 );
            
        for ii=1:rlength
        
            [shape4lin,drshape4lin,dsshape4lin] = iso4lin(rgauss(ii), sgauss(ii));
       
            Jacobiana = jacobiana ( ELX1 ( :, jj ), ELY1 ( : , jj ) , drshape4lin , dsshape4lin );
            invJacob = inv ( Jacobiana );
            detJacob = det ( Jacobiana );
        
            [dxshape4lin,dyshape4lin] = globderiv (invJacob,drshape4lin,dsshape4lin);      
    
            Bmatrica = formiranjeBmatrice_pstrain2 ( dxshape4lin , dyshape4lin );
            Bmatrica2 = formiranjeBmatrice_pstrain2_Eshelby ( dxshape4lin , dyshape4lin );
        
            Klok = Klok  + Bmatrica' * Cmatrica * Bmatrica * wgauss (ii) * detJacob * tgrede ;
            
            shape4lin_save ( (jj-1)*rlength + ii , : ) = shape4lin;
            detJacob_save ( (jj-1)*rlength + ii , : ) = detJacob;
            dxshape4lin_save ( (jj-1)*rlength + ii , : ) = dxshape4lin;
            dyshape4lin_save ( (jj-1)*rlength + ii , : ) = dyshape4lin;
            Bmatrica_save ( (jj-1)*rlength + ii , : ) = reshape(Bmatrica, 1 , 24);
            Bmatrica2_save ( (jj-1)*rlength + ii , : ) = reshape(Bmatrica2, 1 , 32);
                        
        end;
    
    Kl (:,:,jj) = Klok;
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global stiffness matrix calculation

Kg_prije = globmatkrutosti_brzi ( Kl , EL , nC );
Kg = rubni_globmat_brzi ( Kg_prije , Rubni , nC );
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear equation system solver

pomak = Kg \ FV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element strain, stress, Eshelby stress and configurational force
% calculation for Gaussian points

Conf_force_local = zeros ( 8 , nEL );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over elements

for j=1:nEL
    
    Conf_force = zeros ( 8 , 1 );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Local displacements for each element are defined
    
    for k = 1 : 4
        
        Pomaci_elementa ( 2 * k - 1 , 1 ) = pomak ( EL ( 1 + k , j ) * 2 - 1 );
        Pomaci_elementa ( 2 * k , 1 ) = pomak ( EL ( 1 + k , j ) * 2 );
        
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over integration points
    
    for i=1:rlength
        
        hh = ( j - 1 ) * rlength + i;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Local stress and strain vectors with integration point locations

        def_lok = reshape ( Bmatrica_save(hh,:) , 3 , 8 ) * Pomaci_elementa;
        nap_lok = Cmatrica * def_lok;
        
        glob_Gauss_tocke = shape4lin_save ( hh , : ) * [ ELX1(:,j), ELY1(:,j) ];

        Naprezanja ( hh , : ) = [ glob_Gauss_tocke , nap_lok' ];
        Deformacije ( hh , : ) = [ glob_Gauss_tocke , def_lok' ];
%         error('aaa')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deformation gradient and strain energy density in integration
        % point
        
        Def_grad_lista ( hh , : ) = [ 1+dxshape4lin_save(hh,:)*Pomaci_elementa(1:2:8),  dyshape4lin_save(hh,:)*Pomaci_elementa(1:2:8) , dxshape4lin_save(hh,:)*Pomaci_elementa(2:2:8) , 1+dyshape4lin_save(hh,:)*Pomaci_elementa(2:2:8) ];
        W ( hh , : ) = 0.5 * Naprezanja( hh , 3:5 ) * Deformacije( hh , 3:5 )';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eshelby stress in integration point

        Eshelby2 ( hh , 1 ) = W ( hh ) - Naprezanja ( hh , 3 ) * Def_grad_lista ( hh , 1 )  + Naprezanja ( hh , 3 ) - Naprezanja ( hh , 5 ) * Def_grad_lista ( hh , 3 );
        Eshelby2 ( hh , 2 ) = - Naprezanja ( hh , 5 ) * Def_grad_lista ( hh , 1 ) + Naprezanja ( hh , 5 ) - Naprezanja ( hh , 4 ) * Def_grad_lista ( hh , 3 );
        Eshelby2 ( hh , 3 ) = - Naprezanja ( hh , 3 ) * Def_grad_lista ( hh , 2 ) - Naprezanja ( hh , 5 ) * Def_grad_lista ( hh , 4 ) + Naprezanja ( hh , 5 );
        Eshelby2 ( hh , 4 ) = W ( hh ) - Naprezanja ( hh , 5 ) * Def_grad_lista ( hh , 2 ) - Naprezanja ( hh , 4 ) * Def_grad_lista ( hh , 4 ) + Naprezanja ( hh , 4 );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Configuration force Gaussian integration
        
        Conf_force = Conf_force  + reshape(Bmatrica2_save(hh,:),4,8)' * [ Eshelby2((j-1)*rlength+i,1) Eshelby2((j-1)*rlength+i,4) Eshelby2((j-1)*rlength+i,2:3) ]' * wgauss (i) * detJacob_save ( hh ) * tgrede ;


    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Configurational forces for each element - ready for assembling
        
    Conf_force_local (:,j) = Conf_force;
    
end;

% Total potential energy calculation for the whole problem and strain
% energy calculation for each integration point

Energy_FEM = 0.5 * pomak' * Kg * pomak - pomak' * FV;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of configurational forces 

Conf_force = zeros ( 2 * nC , 1 );

for j = 1 : nEL
    
    Conf_force ( 2 * EL ( 2 , j ) - 1 , 1 ) = Conf_force ( 2 * EL ( 2 , j ) - 1 ) + Conf_force_local( 1 , j );
    Conf_force ( 2 * EL ( 2 , j ) , 1 ) = Conf_force ( 2 * EL ( 2 , j ) ) + Conf_force_local( 2 , j );
    Conf_force ( 2 * EL ( 3 , j ) - 1 , 1 ) = Conf_force ( 2 * EL ( 3 , j ) - 1 ) + Conf_force_local( 3 , j );
    Conf_force ( 2 * EL ( 3 , j ) , 1 ) = Conf_force ( 2 * EL ( 3 , j ) ) + Conf_force_local( 4 , j );
    Conf_force ( 2 * EL ( 4 , j ) - 1 , 1 ) = Conf_force ( 2 * EL ( 4 , j ) - 1 ) + Conf_force_local( 5 , j );
    Conf_force ( 2 * EL ( 4 , j ) , 1 ) = Conf_force ( 2 * EL ( 4 , j ) ) + Conf_force_local( 6 , j );
    Conf_force ( 2 * EL ( 5 , j ) - 1 , 1 ) = Conf_force ( 2 * EL ( 5 , j ) - 1 ) + Conf_force_local( 7 , j );
    Conf_force ( 2 * EL ( 5 , j ) , 1 ) = Conf_force ( 2 * EL ( 5 , j ) ) + Conf_force_local( 8 , j );
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing variables

xy_popis { broj_koraka } = xy;
pomak_popis { broj_koraka } = pomak;
Conf_popis { broj_koraka } = Conf_force;
Energy_FEM_popis ( broj_koraka ) = Energy_FEM;
% J_integral_num ( broj_koraka ) = -Conf_force ( Crack_tip * 2 );

if broj_koraka ~= 1
    norma_popis ( broj_koraka ) = norma;
else
    norma = NaN;
end;

if rem ( broj_koraka , 1 ) == 0
    
    disp ( [ 'Broj koraka: ', num2str(broj_koraka) , ' Norma konfiguracijskih sila: ', num2str(norma), ' Energija: ', num2str(Energy_FEM) ] )
    
end;

if broj_koraka == 1 || rem ( broj_koraka , 100 ) == 0
    
    Usporedba_korak = cat(3,Usporedba_korak, broj_koraka);
    Naprezanja_korak = cat(3,Naprezanja_korak , Naprezanja);
    Deformacije_korak = cat(3,Deformacije_korak , Deformacije );
    detJacob_korak = cat(3,detJacob_korak , detJacob_save );
    
end;

end;

% xy_popis = xy_popis ( : , : , 1 : broj_koraka - Oduzmi );
% pomak_popis = pomak_popis ( : , 1 : broj_koraka - Oduzmi );
% Conf_popis = Conf_popis ( : , 1 : broj_koraka - Oduzmi );
% Energy_FEM_popis = Energy_FEM_popis ( : , 1 : broj_koraka - Oduzmi );
% norma_popis = norma_popis ( : , 1 : broj_koraka - Oduzmi );

Naprezanja_korak = cat(2,Naprezanja_korak , zeros(nEL*rlength,4,length(Usporedba_korak)));

Naprezanja_korak ( : , 6 , : ) = (Naprezanja_korak(:,3,:) + Naprezanja_korak(:,4,:))/2 + sqrt( (Naprezanja_korak(:,3,:)-Naprezanja_korak(:,4,:)).^2/4 + Naprezanja_korak(:,5,:).^2 );
Naprezanja_korak ( : , 7 , : ) = (Naprezanja_korak(:,3,:) + Naprezanja_korak(:,4,:))/2 - sqrt( (Naprezanja_korak(:,3,:)-Naprezanja_korak(:,4,:)).^2/4 + Naprezanja_korak(:,5,:).^2 );
Naprezanja_korak ( : , 8 , : ) = atan ( 2 .* Naprezanja_korak(:,5,:) ./ (Naprezanja_korak(:,3,:)-Naprezanja_korak(:,4,:)) ) / 2;
Naprezanja_korak ( : , 9 , : ) = atan ( 2 .* Naprezanja_korak(:,5,:) ./ (Naprezanja_korak(:,3,:)-Naprezanja_korak(:,4,:)) ) / 2 + pi/2;

for ii = 1 : nEL
    
    detJacob_element ( ii , : ) = sum ( detJacob_korak ( (ii - 1 ) * rlength + 1 : ii * rlength ) );
    
end;

% % % % [shape4lin,~,~] = iso4lin(0, 0);
% % % % 
% % % % for ii = 1 : nEL
% % % %     
% % % %     Naprezanja_element ( ii , 1:2 , : )  = shape4lin * [ ELX1(:,ii), ELY1(:,ii) ];
% % % %     Naprezanja_element ( ii , 3 , : ) = sum( Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 6 , : ) .* (detJacob_korak ( (ii - 1) * rlength + 1 : ii * rlength , : ) ) / detJacob_element ( ii ) );
% % % %     Naprezanja_element ( ii , 4 , : ) = sum( Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 7 , : ) .* (detJacob_korak ( (ii - 1) * rlength + 1 : ii * rlength , : ) ) / detJacob_element ( ii ) );
% % % %     Naprezanja_element ( ii , 5 , : ) = sum( Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 8 , : ) .* (detJacob_korak ( (ii - 1) * rlength + 1 : ii * rlength , : ) ) / detJacob_element ( ii ) );
% % % %     Naprezanja_element ( ii , 6 , : ) = sum( Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 9 , : ) .* (detJacob_korak ( (ii - 1) * rlength + 1 : ii * rlength , : ) ) / detJacob_element ( ii ) );
% % % % 
% % % % end;


% Stress recovery procedure from Boulder Colorado IFEM chapter 28
% Stress interpolation from Gauss points to nodal points with additional
% nodal stress smoothening

Stress_recovery = zeros ( nEL*rlength , 5 , length(Usporedba_korak) );
Recovery_matrix = [ 1+sqrt(3)/2 -1/2 1-sqrt(3)/2 -1/2; -1/2 1+sqrt(3)/2 -1/2 1-sqrt(3)/2; 1-sqrt(3)/2 -1/2 1+sqrt(3)/2 -1/2; -1/2 1-sqrt(3)/2 -1/2 1+sqrt(3)/2 ];

for jj = 1 : length(Usporedba_korak)

for ii = 1 : nEL
    
    Stress_recovery ( (ii-1)*rlength+1:ii*rlength , 1 , jj ) = Recovery_matrix * Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 3 , jj );
    Stress_recovery ( (ii-1)*rlength+1:ii*rlength , 2 , jj ) = Recovery_matrix * Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 4 , jj );
    Stress_recovery ( (ii-1)*rlength+1:ii*rlength , 3 , jj ) = Recovery_matrix * Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 5 , jj );
    Stress_recovery ( (ii-1)*rlength+1:ii*rlength , 4 , jj ) = Recovery_matrix * Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 6 , jj );
    Stress_recovery ( (ii-1)*rlength+1:ii*rlength , 5 , jj ) = Recovery_matrix * Naprezanja_korak ( (ii-1)*rlength+1:ii*rlength , 7 , jj );

end;

end;

Stress_point = zeros ( nC , 5 , length(Usporedba_korak) );

for jj = 1 : length(Usporedba_korak)

for ii = 1 : nEL
    
    Stress_point ( EL ( 2 , ii ) , 1:5 , jj ) = Stress_point ( EL ( 2 , ii ) , 1:5 , jj ) + Stress_recovery ( (ii-1)*rlength+1 , : , jj );
    Stress_point ( EL ( 3 , ii ) , 1:5 , jj ) = Stress_point ( EL ( 3 , ii ) , 1:5 , jj ) + Stress_recovery ( (ii-1)*rlength+2 , : , jj );
    Stress_point ( EL ( 4 , ii ) , 1:5 , jj ) = Stress_point ( EL ( 4 , ii ) , 1:5 , jj ) + Stress_recovery ( (ii-1)*rlength+3 , : , jj );
    Stress_point ( EL ( 5 , ii ) , 1:5 , jj ) = Stress_point ( EL ( 5 , ii ) , 1:5 , jj ) + Stress_recovery ( (ii-1)*rlength+4 , : , jj );
    
end;

end;

for ii = 1 : nC
    
    Node_occurence ( ii , 1 ) = sum ( sum ( ii == EL ( 2:5 , : ) ) );
    
end;

for ii = 1 : nC
    
    Stress_point ( ii , : , : ) = Stress_point ( ii , : , : ) ./ Node_occurence ( ii );

end;


ELX(2:5,:) = ELX1;
ELY(2:5,:) = ELY1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some fracture mechanics characteristics

% Analytical results

% Stress_intensity_factor_an = ( 1.12 - 0.23 * (notch/hgrede) + 10.56 * (notch/hgrede)^2 - 21.74 * (notch/hgrede)^3 + 30.42 * (notch/hgrede)^4 ) * p0 * sqrt( pi * notch );
% Stress_intensity_factor_an2 = ( 1.12 - 0.231 * (notch/hgrede) + 10.55 * (notch/hgrede)^2 - 21.72 * (notch/hgrede)^3 + 30.39 * (notch/hgrede)^4 ) * p0 * sqrt( pi * notch );
Stress_intensity_factor_an3 = ( 1.12 - 0.23 * (notch/hgrede) + 10.6 * (notch/hgrede)^2 - 21.7 * (notch/hgrede)^3 + 30.4 * (notch/hgrede)^4 ) * p0 * sqrt( pi * notch );
% 
% J_integral_an = Stress_intensity_factor_an^2 / E * ( 1 - Nu^2 );
% J_integral_an2 = Stress_intensity_factor_an2^2 / E * ( 1 - Nu^2 );
J_integral_an3 = Stress_intensity_factor_an3^2 * ( 1 - Nu ) / ( 2 * G );
% 
% 
% Normalised_J_an = J_integral_an * Lame_1 / ( p0^2 * notch );
% Normalised_J_an2 = J_integral_an2 * Lame_1 / ( p0^2 * notch );
% Normalised_J_an3 = J_integral_an3 * Lame_1 / ( p0^2 * notch );
% Normalised_J_num = J_integral_num * Lame_1 / ( p0^2 * notch );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results visualisation - needs improvement

[ xyDEF , xyPLOT ] = uredjenje_deformacija ( xy , pomak );

for i = 1:nEL
    ELXPLOT2 ( 1 , i ) = EL ( 1 , i );
    ELXPLOT2 ( 2 , i ) = xyPLOT ( 1, EL ( 2 , i ) );
    ELXPLOT2 ( 3 , i ) = xyPLOT ( 1, EL ( 3 , i ) );
    ELXPLOT2 ( 5 , i ) = xyPLOT ( 1, EL ( 5 , i ) );
    ELXPLOT2 ( 4 , i ) = xyPLOT ( 1, EL ( 4 , i ) );
end;
for i = 1:nEL
    ELYPLOT2 ( 1 , i ) = EL ( 1 , i );
    ELYPLOT2 ( 2 , i ) = xyPLOT ( 2, EL ( 2 , i ) );
    ELYPLOT2 ( 3 , i ) = xyPLOT ( 2, EL ( 3 , i ) );
    ELYPLOT2 ( 5 , i ) = xyPLOT ( 2, EL ( 5 , i ) );
    ELYPLOT2 ( 4 , i ) = xyPLOT ( 2, EL ( 4 , i ) );
end;

ELX1PLOT = ELXPLOT2 ( 2:5 , : );
ELY1PLOT = ELYPLOT2 ( 2:5 , : );

% figure
% hold on
% plot(ELXPLOT,ELYPLOT,'k-')

% [ELX , ELY , EL ] = formiranje_4kutnihelemenata( xy, nlelemenata, nhelemenata );

ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT1 = [ELX(2:5,:);ELX(2,:)];
ELYPLOT1 = [ELY(2:5,:);ELY(2,:)];

% plot(ELXPLOT1,ELYPLOT1,'r--')
% hold off

figure 
hold on

for i = 1:nEL
    ELXPLOT3 ( 1 , i ) = xyDEF ( 1, EL ( 2 , i ) );
    ELXPLOT3 ( 2 , i ) = xyDEF ( 1, EL ( 3 , i ) );
    ELXPLOT3 ( 4 , i ) = xyDEF ( 1, EL ( 5 , i ) );
    ELXPLOT3 ( 3 , i ) = xyDEF ( 1, EL ( 4 , i ) );
    ELXPLOT3 ( 5 , i ) = xyDEF ( 1, EL ( 2 , i ) );
end;
for i = 1:nEL
    ELYPLOT3 ( 1 , i ) = xyDEF ( 2, EL ( 2 , i ) );
    ELYPLOT3 ( 2 , i ) = xyDEF ( 2, EL ( 3 , i ) );
    ELYPLOT3 ( 4 , i ) = xyDEF ( 2, EL ( 5 , i ) );
    ELYPLOT3 ( 3 , i ) = xyDEF ( 2, EL ( 4 , i ) );
    ELYPLOT3 ( 5 , i ) = xyDEF ( 2, EL ( 2 , i ) );
end;

plot(ELXPLOT3,ELYPLOT3,'k-')
% plot(ELXPLOT1,ELYPLOT1,'r--')
hold off

% Q4_Conf_force_to_VTK(nC,nEL,xy_popis(:,:,reshape(Usporedba_korak,1,length(Usporedba_korak))),EL,pomak_popis(:,reshape(Usporedba_korak,1,length(Usporedba_korak))),length(Usporedba_korak),Stress_point)












