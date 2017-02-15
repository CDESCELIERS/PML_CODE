clear
clc
close all

load('num_ddl_u.ascii');
num_ddl_u

load('num_ddl_E.ascii');
num_ddl_E

disp('Table_Elements');
load('Table_Elements.ascii');
Ne = size(Table_Elements,1);

disp('Table_Noeuds');
load('Table_Noeuds.ascii');
Nn = size(Table_Noeuds,2);

disp('M_DDL_U');
load('DATA_M_DDL_U.ascii');
M_DDL_U = sparse(DATA_M_DDL_U(:,2), DATA_M_DDL_U(:,3), DATA_M_DDL_U(:,1), Nn*2, num_DDL_u);

ddl = 1:2*Nn;
ddl_fixed = ddl( M_DDL_U * ones(num_DDL_u,1) == 0 )';

disp('M_DDL_E');
load('DATA_M_DDL_E.ascii');
M_DDL_E = sparse(DATA_M_DDL_E(:,2), DATA_M_DDL_E(:,3), DATA_M_DDL_E(:,1), Nn*3, num_DDL_E);
ddl = 1:3*Nn;
Null_E = ddl( M_DDL_E * ones(num_DDL_E,1) == 0 )';

IDDL_U = 1:(2*Nn);
IDDL_U(ddl_fixed)= [];

IDDL_E = 1:3*Nn;
IDDL_E(Null_E) = [];

Ustock_pml = zeros(2*Nn,1);
nddl = num_DDL_u;
nddl_pml = num_DDL_E;

nddl_glob =  num_DDL_u + num_DDL_E;

disp('DATA_M_K_dyn');
load('DATA_M_K_dyn.ascii');
M_K_dyn = sparse(DATA_M_K_dyn(:,3), DATA_M_K_dyn(:,4), DATA_M_K_dyn(:,1)+1i* DATA_M_K_dyn(:,2), nddl_glob,nddl_glob);

load('DATA_V_F_pml_complex.ascii');
V_F_pml =  DATA_V_F_pml_complex(:,1) + 1i*DATA_V_F_pml_complex(:,2) ;

tic
U_pml = M_K_dyn\V_F_pml;

Ustock_pml(IDDL_U) = U_pml(1:nddl);
toc

%%%%%%  VTK FILE FOR THE RESULTS %%%%%%

filename = 'Resultats_freq_pml';

disp(['Ecriture du fichiers vtk ' filename]);

s = sprintf ('%s.vtk', filename);
fid = fopen(s,'w');

if (fid<0)
    error('CANNOT OPEN OUTPUT FILE');
end

fprintf ( fid, '# vtk DataFile Version 2.0\n' );

fprintf ( fid, 'Resultats etude PML en frequence \n' );

fprintf ( fid, 'ASCII\n' );
fprintf ( fid, '\n' );
fprintf ( fid, 'DATASET UNSTRUCTURED_GRID\n' );
fprintf ( fid, 'POINTS %d double\n', Nn);
for k=1:Nn
    fprintf ( fid, '%f %f %f\n', Table_Noeuds(1,k), Table_Noeuds(2,k), 0);
end
fprintf ( fid, '\n' );
fprintf ( fid, 'CELLS %d %d\n', Ne, 7*Ne);
for e=1:Ne
    bid = Table_Elements(e,2:end);
    fprintf ( fid, '6 %d %d %d %d %d %d\n', bid-1);
end
fprintf ( fid, '\n' );
fprintf ( fid, 'CELL_TYPES %d\n',Ne);
for e=1:Ne
    fprintf ( fid, '22\n');
end
fprintf ( fid, '\n' );
fprintf ( fid, 'POINT_DATA %d\n', Nn);

fprintf ( fid, 'VECTORS Disp double');
fprintf ( fid, '    %e  %e  0 \n',Ustock_pml(:));
fprintf ( fid, '\n' );

fclose(fid);












