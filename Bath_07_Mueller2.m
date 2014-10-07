% Configurational forces 2d problem - Bath 06.10.2014
% Mesh improvement by configurational force analysis. Example problem from
% Mueller R., Kolling S., Gross D., "On configurational forces in the
% context of the finite element method" - example 4.1.2.Clamped block under
% tension. Procedure author: Nikola Lustig, mag.ing.aedif.

clear; 
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters

lgrede = 10;
hgrede = 10;
tgrede = 1;

nlelemenata = 6;
nhelemenata = 6;

E = 5*1000;
Nu = 0.25;

% Lame constants calculation for comparison with paper
G = E/(2*(1+Nu));
Lame_1 = E * Nu / ( ( 1 + Nu ) * ( 1 - 2 * Nu ) );

% p0 defines a distributed load
p0 = 1000;
% Sila defines a concentrated load
% Sila = -100;

% Tangential movement of boundary nodes is allowed everywhere apart from
% the loaded edge
tang_remesh = 1;

c = 0.01;
max_koraka = 10000;
tol = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite element mesh creation

xy = pravokutni_mesh4kutni( lgrede , hgrede , nlelemenata , nhelemenata );
[ELX , ELY , EL ] = formiranje_4kutnihelemenata( xy, nlelemenata, nhelemenata );

ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT = [ELX(2:5,:);ELX(2,:)];
ELYPLOT = [ELY(2:5,:);ELY(2,:)];


nC = length ( xy );
nEL = nlelemenata * nhelemenata;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining non-boundary nodes which are considered in configurational
% force based mesh improvement

Internal_node = [];
External_node = [];
Tangential_node_x = [];
Tangential_node_y = [];
Fixed_node = [];

for ii = 1 : nC
    
    if xy ( 1 , ii ) ~= 0 && xy ( 1 , ii ) ~= lgrede && xy ( 2 , ii ) ~= 0 && xy ( 2 , ii ) ~= hgrede
        
        Internal_node = [ Internal_node ii ];
        
    end;
    
end;

if tang_remesh == 1
    
    for ii = 1 : nC
        
        if xy ( 1 , ii ) == 0  
            
            Tangential_node_x = [ Tangential_node_x ii ];
            
        end;
        
        if ( xy ( 2 , ii ) == 0 || xy ( 2 , ii ) == hgrede ) && xy ( 1 , ii ) ~= lgrede
            
            Tangential_node_y = [ Tangential_node_y ii ];
            
        end;
        
    end;
    
end;
    
for ii = 1 : nC
    
    if ~(any(ii==Internal_node) || any(ii==Tangential_node_x) || any(ii==Tangential_node_y))
        
        Fixed_node = [ Fixed_node ii ];
        
    end;
    
end;
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plane strain state for isotropic linear elastic material
Smatrica = [ (1-Nu^2)/E -Nu*(1+Nu)/E 0; -Nu*(1+Nu)/E (1-Nu^2)/E 0; 0 0 1/G ];
Cmatrica = inv(Smatrica);

% % Plane stress state for isotropic linear elastic material
% Cmatrica = ( E / ( 1 - Nu^2 ) ) * [ 1 Nu 0; Nu 1 0; 0 0 (1-Nu)/2 ]; 


% Gaussian numerical quadrature rule - 2x2 integration
[ rgauss, sgauss, wgauss ] = gaussquad2_shift;
% [ rgauss, sgauss, wgauss ] = gaussquad2;
rlength = length ( rgauss );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary values and loading definition

Rubni = 1 : nlelemenata+1 : nC;
Rubni = [ Rubni ; zeros(1,length(Rubni)) ; zeros(1,length(Rubni)) ];


FV = zeros ( 2 * nC , 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a vertical end-node concentrated load
% FV ( 2 * nC ) = Sila;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a vertical distributed load
% Sila = p0 / (nlelemenata+1);
% FV ( 2*(nC - (1:nlelemenata+1)) ) = Sila;
% FV(2*nC) = Sila/2;
% FV(2*(nC-nlelemenata))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this for a horizontal distributed load
Sila = p0 / (nhelemenata);
FV ( 2*(nlelemenata+1:nlelemenata+1:nC) - 1 ) = Sila;
FV(2*nC - 1) = Sila/2;
FV(2*(nlelemenata+1)-1)=Sila/2;
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
xy_popis = zeros ( 2 , nC , max_koraka );
pomak_popis = zeros ( nC * 2 , max_koraka );
Conf_popis = zeros ( nC * 2 , max_koraka );
Energy_FEM_popis = zeros ( 1 , max_koraka );
norma_popis = zeros ( 1 , max_koraka );
Oduzmi = 0;

Usporedba_korak = [  ];
Naprezanja_korak = [  ];
Deformacije_korak = [  ];
detJacob_korak = [  ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of loop over elements for mesh improvement - NODE SHIFTING

for broj_koraka = 1 : max_koraka
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Checking convergence criteria and mesh adaptation
    
    if broj_koraka ~= 1
        
        Cforce_operational = Conf_force;
        
        Cforce_operational ( [2*Tangential_node_x-1,2*Tangential_node_y,2*Fixed_node,2*Fixed_node-1], 1 ) = 0;    
            
        norma = norm ( Cforce_operational );
        
        if norma <= tol
            
            disp('Konvergencija')
            Oduzmi = 1;

    
            Usporedba_korak = cat(3,Usporedba_korak, broj_koraka-1);
            Naprezanja_korak = cat(3,Naprezanja_korak , Naprezanja);
            Deformacije_korak = cat(3,Deformacije_korak , Deformacije );
            detJacob_korak = cat(3,detJacob_korak , detJacob_save );
    

            break
            
        end;
        
        xy  = xy  - c * reshape ( Cforce_operational , 2 , length ( xy ) );
        
        for i = 1:nEL
            ELX1 ( 1 , i ) = xy ( 1, EL ( 2 , i ) );
            ELX1 ( 2 , i ) = xy ( 1, EL ( 3 , i ) );
            ELX1 ( 4 , i ) = xy ( 1, EL ( 5 , i ) );
            ELX1 ( 3 , i ) = xy ( 1, EL ( 4 , i ) );
            ELY1 ( 1 , i ) = xy ( 2, EL ( 2 , i ) );
            ELY1 ( 2 , i ) = xy ( 2, EL ( 3 , i ) );
            ELY1 ( 4 , i ) = xy ( 2, EL ( 5 , i ) );
            ELY1 ( 3 , i ) = xy ( 2, EL ( 4 , i ) );
        end;
        
    end;
    
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

% Calculation of the Eshelby stress in each integration point

% Approach from Mueller2002 - linear elastic - take one
% % % Eshelby = zeros ( rlength * nEL , 1 );
% % % 
% % % for ii = 1 : rlength * nEL
% % %     
% % %     Eshelby ( ii , 1 ) = W ( ii ) - Naprezanja ( ii , 3 ) * ( Def_grad_lista ( ii , 1 ) - 1 ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 3 );
% % %     Eshelby ( ii , 2 ) = - Naprezanja ( ii , 3 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 5 ) * ( Def_grad_lista ( ii , 4 ) - 1 );
% % %     Eshelby ( ii , 3 ) = - Naprezanja ( ii , 5 ) * ( Def_grad_lista ( ii , 1 ) - 1 ) - Naprezanja ( ii , 4 ) * Def_grad_lista ( ii , 3 );
% % %     Eshelby ( ii , 4 ) = W ( ii ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 4 ) * ( Def_grad_lista ( ii , 4 ) - 1 );
% % % 
% % % end;

% Approach from Mueller2002 - linear elastic - take two

% % % Eshelby2 = zeros ( rlength * nEL , 1 );
% % % 
% % % for ii = 1 : rlength * nEL
% % %     
% % %     Eshelby2 ( ii , 1 ) = W ( ii ) - Naprezanja ( ii , 3 ) * Def_grad_lista ( ii , 1 )  + Naprezanja ( ii , 3 ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 3 );
% % %     Eshelby2 ( ii , 2 ) = - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 1 ) + Naprezanja ( ii , 5 ) - Naprezanja ( ii , 4 ) * Def_grad_lista ( ii , 3 );
% % %     Eshelby2 ( ii , 3 ) = - Naprezanja ( ii , 3 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 4 ) + Naprezanja ( ii , 5 );
% % %     Eshelby2 ( ii , 4 ) = W ( ii ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 4 ) * Def_grad_lista ( ii , 4 ) + Naprezanja ( ii , 4 );
% % % 
% % % end;

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

xy_popis ( : , : , broj_koraka ) = xy;
pomak_popis ( : , broj_koraka ) = pomak;
Conf_popis ( : , broj_koraka ) = Conf_force;
Energy_FEM_popis ( broj_koraka ) = Energy_FEM;

if broj_koraka ~= 1
    norma_popis ( broj_koraka ) = norma;
end;

if rem ( broj_koraka , 20 ) == 0
    
    disp ( [ 'Broj koraka: ', num2str(broj_koraka) , ' Norma konfiguracijskih sila: ', num2str(norma), ' Energija: ', num2str(Energy_FEM) ] )
    
end;

if broj_koraka == 1 || rem ( broj_koraka , 100 ) == 0
    
    Usporedba_korak = cat(3,Usporedba_korak, broj_koraka);
    Naprezanja_korak = cat(3,Naprezanja_korak , Naprezanja);
    Deformacije_korak = cat(3,Deformacije_korak , Deformacije );
    detJacob_korak = cat(3,detJacob_korak , detJacob_save );
    
end;

end;

xy_popis = xy_popis ( : , : , 1 : broj_koraka - Oduzmi );
pomak_popis = pomak_popis ( : , 1 : broj_koraka - Oduzmi );
Conf_popis = Conf_popis ( : , 1 : broj_koraka - Oduzmi );
Energy_FEM_popis = Energy_FEM_popis ( : , 1 : broj_koraka - Oduzmi );
norma_popis = norma_popis ( : , 1 : broj_koraka - Oduzmi );

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

figure
hold on
plot(ELXPLOT,ELYPLOT,'k-')

[ELX , ELY , EL ] = formiranje_4kutnihelemenata( xy, nlelemenata, nhelemenata );

ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT1 = [ELX(2:5,:);ELX(2,:)];
ELYPLOT1 = [ELY(2:5,:);ELY(2,:)];

plot(ELXPLOT1,ELYPLOT1,'r--')
hold off

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
plot(ELXPLOT1,ELYPLOT1,'r--')
hold off

Q4_Conf_force_to_VTK(nC,nEL,xy_popis(:,:,reshape(Usporedba_korak,1,length(Usporedba_korak))),EL,pomak_popis(:,reshape(Usporedba_korak,1,length(Usporedba_korak))),length(Usporedba_korak),Stress_point)












