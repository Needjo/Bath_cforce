% Configurational forces 2d problem - Bath 22.09.2014

clear; 
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters

lgrede = 1;
hgrede = 1;
tgrede = 1;

nlelemenata = 1;
nhelemenata = 1;

E = 1000;
Nu = 0.3;

G = E/(2*(1+Nu));

Sila = 200;

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
% Local stiffness matrix formulation

% Plane strain state for istoropic linear elastic material
Smatrica = [ (1-Nu^2)/E -Nu*(1+Nu)/E 0; -Nu*(1+Nu)/E (1-Nu^2)/E 0; 0 0 1/G ];
Cmatrica = inv(Smatrica);

% % Plane stress state for istoropic linear elastic material
% Cmatrica = ( E / ( 1 - Nu^2 ) ) * [ 1 Nu 0; Nu 1 0; 0 0 (1-Nu)/2 ]; 


% Gaussian numerical quadrature rule - 2x2 integration
[ rgauss, sgauss, wgauss ] = gaussquad2;
rlength = length ( rgauss );


% Local stiffness matrix calculation - linear elastic isotropic elements
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
        
            Klok = Klok  + Bmatrica' * Cmatrica * Bmatrica * wgauss (ii) * detJacob * tgrede ;
                        
        end;
    
    % Postavljanje gotove lokalne matrice krutosti na redno mjesto
    % elementa. Dobivamo listu lokalnih matrica krutosti dugacku kao ukupni
    % broj konacnih elemenata
    
    Kl (:,:,jj) = Klok;
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary values and loading definition
% Rubni = rubniuvjet_ptest ( xy , nC );

% Rubni = [1   ,  2 ,    3   ,  4;
%      0  ,   0    , 0    , 0;
%      0   ,  0   ,  1   ,  1];

% Rubni = [ 1 2 3; 0 1 1 ; 0 0 0]

Rubni = [];
for ii = 1 : nC
    
    if xy ( 1 , ii ) == 0
        
        Rubni = [ Rubni ii ];
        
    end;
    
end;
  
Rubni = [ Rubni ; zeros(1,length(Rubni)) ; ones(1,length(Rubni)) ];
Rubni ( 3 , 1 ) = 0;


% Rubni = [ 1 2 3 4; 0 1 0 1 ; 0 0 0 0];

% Rubni = [ 1 2 3 4; 0 0 0 0 ; 0 0 1 1];

% Rubni = [ 1 2 3 4; 0 1 0 1 ; 0 0 1 1];

Rubni = [ 1 2 3 4; 0 0 1 1 ; 0 0 1 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global stiffness matrix calculation

Kg_prije = globmatkrutosti_brzi ( Kl , EL , nC );
Kg = rubni_globmat_brzi ( Kg_prije , Rubni , nC );

FV = zeros ( 2 * nC , 1 );
for ii = 1 : nC
    if xy ( 2 , ii ) == hgrede && ( xy ( 1 , ii ) == 0 || round(xy ( 1 , ii )*1000)/1000 == lgrede )
        FV ( 2 * ii , 1 ) = Sila / (2*( nlelemenata ));
    elseif xy ( 2 , ii ) == hgrede
        FV ( 2 * ii , 1 ) = Sila / ( nlelemenata );
    end;
end; 

% FV = zeros ( 2 * nC , 1 );
% FV (end - 4 ) = Sila/2;
% FV(end - 2) = Sila;
% FV ( end ) = Sila/2;

FV = zeros ( 2 * nC , 1 );
for ii = 1 : nC
    if round(xy ( 1 , ii )*1000)/1000 == lgrede && ( xy ( 2 , ii ) == 0 || round(xy ( 2 , ii )*1000)/1000 == hgrede )
        FV ( 2 * ii-1 , 1 ) = Sila / (2*( nhelemenata ));
    elseif round(xy ( 1 , ii )*1000)/1000 == lgrede
        FV ( 2 * ii-1 , 1 ) = Sila / ( nhelemenata );
    end;
end;

% FV = [ 0 0 0 0 0 Sila/2 0 Sila/2 ]';

FV = [ 0 0 0 0 0 Sila/2 0 Sila/2 ]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear equation system solver

pomak = Kg \ FV;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element strain and stress calculation for Gaussian points



% Petlja za pozivanje svakog elementa
napx_lista=[];
napy_lista=[];
napxy_lista=[];
defx_lista=[];
defy_lista=[];
defxy_lista=[];
Gauss_tocke=[];

for j=1:nEL
    Pomaci_elementa = [];
    Pomaci_elementa_x = [];
    Pomaci_elementa_y = [];
    for k=1:4
        Pomaci_elementa = [ Pomaci_elementa; pomak(EL(1+k,j)*2-1) ];
        Pomaci_elementa = [ Pomaci_elementa; pomak(EL(1+k,j)*2) ];
        Pomaci_elementa_x = [Pomaci_elementa_x;pomak(EL(1+k,j)*2-1) ];
        Pomaci_elementa_y = [Pomaci_elementa_y;pomak(EL(1+k,j)*2) ];
    end;

      Def_grad_lista = [];
      
    % Petlja za provoðenje Gauss-ove integracije
    
    for i=1:rlength
         
        % Procedura koja poziva funckije oblika i njihove derivacije za
        % pojedinu vrstu konaènih elemenata
        
        [shape4lin,drshape4lin,dsshape4lin] = iso4lin(rgauss(i), sgauss(i));
        
        % Procedura koja izraèunava Jacobi-evu matricu za transformaciju iz
        % lokalnih derivacija u globalne
        
        Jacobiana = jacobiana ( ELX1 ( :, j ), ELY1 ( : , j ) , drshape4lin , dsshape4lin );
        invJacob = inv ( Jacobiana );
        detJacob = det ( Jacobiana );
        
        % Procedura koja izraèunava globalne derivacije pomoæu Jacobi-eve
        % matrice i lokalnih derivacija funkcija oblika
        
        [dxshape4lin,dyshape4lin]=globderiv(invJacob,drshape4lin,dsshape4lin);
        
        % Procedura koja odreðuje B matricu sustava - matricu veze
        % deformacija i pomaka
        
        Bmatrica = formiranjeBmatrice_pstrain2 ( dxshape4lin , dyshape4lin );
        
        def_lok = Bmatrica * Pomaci_elementa;
        nap_lok = Cmatrica * def_lok;
        
        glob_Gauss_tocke = shape4lin * [ ELX(2:5,j), ELY(2:5,j) ];
        
        napx_lista = [ napx_lista; nap_lok(1) ];
        napy_lista = [ napy_lista; nap_lok(2) ];
        napxy_lista = [ napxy_lista; nap_lok(3) ];
        defx_lista = [ defx_lista; def_lok(1) ];
        defy_lista = [ defy_lista; def_lok(2) ];
        defxy_lista = [ defxy_lista; def_lok(3) ];
        Gauss_tocke = [ Gauss_tocke; glob_Gauss_tocke ];
        
        Def_grad_lista = [ Def_grad_lista ; 1+dxshape4lin*Pomaci_elementa_x,  dyshape4lin*Pomaci_elementa_x,dxshape4lin*Pomaci_elementa_y, 1+dyshape4lin*Pomaci_elementa_y  ];
        
    end;
end;
Naprezanja = [ Gauss_tocke, napx_lista, napy_lista, napxy_lista ];
Deformacije = [ Gauss_tocke, defx_lista, defy_lista, defxy_lista ];

for ii = 1 : rlength * nEL
    
    Naprezanja_tenzor ( : , : , ii ) = [Naprezanja(ii,3),Naprezanja(ii,5);Naprezanja(ii,5),Naprezanja(ii,4)];
    Deformacije_tenzor ( : , : , ii ) = [Deformacije(ii,3),Deformacije(ii,5);Deformacije(ii,5),Deformacije(ii,4)];
    Deformacijski_gradijent ( : , : , ii ) = [Def_grad_lista(ii,1),Def_grad_lista(ii,2);Def_grad_lista(ii,3),Def_grad_lista(ii,4)];
    J ( ii ) = det ( Deformacijski_gradijent ( : , : , ii ) );
        
end;


% [ Naprezanja , Deformacije ] = naprezanja_deformacije ( pomak , EL , ELX , ELY , nEL , Cmatrica );





% Total potential energy calculation for the whole problem and strain
% energy calculation for each integration point
Energy_FEM = 0.5 * pomak' * Kg * pomak - pomak' * FV;

W = zeros ( rlength * nEL , 1 );
for ii = 1 : rlength * nEL
    
    W ( ii ) = 0.5 * [Naprezanja(ii,3:5)] * [Deformacije(ii,3:5)]';
% ,Naprezanja(ii,5)         ,Deformacije(ii,5)
end;

% Calculation of the Eshelby stress in each integration point

% Approach directly from 1d - very likely to be wrong
Eshelby_tenzor = zeros ( 2 , 2 , rlength * nEL );
for ii = 1 : rlength * nEL
    
    Eshelby_tenzor ( : , : , ii ) = W ( ii ) .* eye ( 2 ) -  Naprezanja_tenzor ( : , : , ii ) * Deformacije_tenzor ( : , : , ii );
    
end;

% % Approach from Mueller2002 - hyperelastic - rejected due to using a Cauchy
% % stress in the first place to define the strain energy density
% Eshelby_tenzor2 = zeros ( 2 , 2 , rlength * nEL );
% for ii = 1 : rlength * nEL
%     
%     Eshelby_tenzor2 ( : , : , ii ) = W ( ii ) .* eye ( 2 ) -  J ( ii ) * Deformacijski_gradijent ( : , : , ii )' * Naprezanja_tenzor ( : , : , ii ) * inv ( Deformacijski_gradijent ( : , : , ii )' );
%     
% end;

% Approach from Mueller2002 - hyperelastic - not using all the
% transormations, just the basic facts, can you show me where it hurts,
% there is no pain you are receding, no distant ship smoke on the horizon
Eshelby_tenzor3 = zeros ( 2 , 2 , rlength * nEL );
for ii = 1 : rlength * nEL
    
    Eshelby_tenzor3 ( : , : , ii ) = W ( ii ) .* eye ( 2 ) -  Naprezanja_tenzor ( : , : , ii ) * ( Deformacijski_gradijent ( : , : , ii ) - eye ( 2 ) );
    
end;


% Approach from Mueller2002 - linear elastic - take one
Eshelby = zeros ( rlength * nEL , 1 );

for ii = 1 : rlength * nEL
    
    Eshelby ( ii , 1 ) = W ( ii ) - Naprezanja ( ii , 3 ) * ( Def_grad_lista ( ii , 1 ) - 1 ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 3 );
    Eshelby ( ii , 2 ) = - Naprezanja ( ii , 3 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 5 ) * ( Def_grad_lista ( ii , 4 ) - 1 );
    Eshelby ( ii , 3 ) = - Naprezanja ( ii , 5 ) * ( Def_grad_lista ( ii , 1 ) - 1 ) - Naprezanja ( ii , 4 ) * Def_grad_lista ( ii , 3 );
    Eshelby ( ii , 4 ) = W ( ii ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 4 ) * ( Def_grad_lista ( ii , 4 ) - 1 );

end;

% Approach from Mueller2002 - linear elastic - take two

Eshelby2 = zeros ( rlength * nEL , 1 );

for ii = 1 : rlength * nEL
    
    Eshelby2 ( ii , 1 ) = W ( ii ) - Naprezanja ( ii , 3 ) * Def_grad_lista ( ii , 1 )  + Naprezanja ( ii , 3 ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 3 );
    Eshelby2 ( ii , 2 ) = - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 1 ) + Naprezanja ( ii , 5 ) - Naprezanja ( ii , 4 ) * Def_grad_lista ( ii , 3 );
    Eshelby2 ( ii , 3 ) = - Naprezanja ( ii , 3 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 4 ) + Naprezanja ( ii , 5 );
    Eshelby2 ( ii , 4 ) = W ( ii ) - Naprezanja ( ii , 5 ) * Def_grad_lista ( ii , 2 ) - Naprezanja ( ii , 4 ) * Def_grad_lista ( ii , 4 ) + Naprezanja ( ii , 4 );

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Result visualisation - needs improvement

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

ELXPLOT1 = [ ELX1PLOT; ELXPLOT2(2,:)];
ELYPLOT1 = [ ELY1PLOT; ELYPLOT2(2,:)];
plot(ELXPLOT1,ELYPLOT1,'r--')
hold off
















