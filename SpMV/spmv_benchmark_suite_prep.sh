#!/bin/bash

cd utilities
make
cd ..
mkdir matrices
cd matrices

echo "Fetching benchmark matrices"

wget https://sparse.tamu.edu/MM/Bourchtein/atmosmodd.tar.gz
wget https://sparse.tamu.edu/MM/Wissgott/parabolic_fem.tar.gz
wget https://sparse.tamu.edu/MM/Rajat/rajat30.tar.gz
wget https://sparse.tamu.edu/MM/Bodendiek/CurlCurl_3.tar.gz
wget https://sparse.tamu.edu/MM/Um/offshore.tar.gz
wget https://sparse.tamu.edu/MM/Botonakis/FEM_3D_thermal2.tar.gz
wget https://sparse.tamu.edu/MM/Schenk/nlpkkt80.tar.gz
wget https://sparse.tamu.edu/MM/PARSEC/CO.tar.gz
wget https://sparse.tamu.edu/MM/Dziekonski/gsm_106857.tar.gz
wget https://sparse.tamu.edu/MM/INPRO/msdoor.tar.gz
wget https://sparse.tamu.edu/MM/GHS_indef/bmw3_2.tar.gz
wget https://sparse.tamu.edu/MM/BenElechi/BenElechi1.tar.gz
wget https://sparse.tamu.edu/MM/Oberwolfach/t3dh.tar.gz
wget https://sparse.tamu.edu/MM/Koutsovasilis/F2.tar.gz
wget https://sparse.tamu.edu/MM/Williams/consph.tar.gz
wget https://sparse.tamu.edu/MM/PARSEC/SiO2.tar.gz
wget https://sparse.tamu.edu/MM/Norris/torso1.tar.gz
wget https://sparse.tamu.edu/MM/Dziekonski/dielFilterV3real.tar.gz
wget https://sparse.tamu.edu/MM/Fluorem/RM07R.tar.gz
wget https://sparse.tamu.edu/MM/DNVS/m_t1.tar.gz
wget https://sparse.tamu.edu/MM/GHS_psdef/crankseg_2.tar.gz
wget https://sparse.tamu.edu/MM/ND/nd24k.tar.gz
wget https://sparse.tamu.edu/MM/TSOPF/TSOPF_RS_b2383.tar.gz
wget https://sparse.tamu.edu/MM/Belcastro/mouse_gene.tar.gz
wget https://sparse.tamu.edu/MM/Belcastro/human_gene1.tar.gz

echo "Unpacking matrices"

tar -xvf atmosmodd.tar.gz
tar -xvf parabolic_fem.tar.gz
tar -xvf rajat30.tar.gz
tar -xvf CurlCurl_3.tar.gz
tar -xvf offshore.tar.gz
tar -xvf FEM_3D_thermal2.tar.gz
tar -xvf nlpkkt80.tar.gz
tar -xvf CO.tar.gz
tar -xvf gsm_106857.tar.gz
tar -xvf msdoor.tar.gz
tar -xvf bmw3_2.tar.gz
tar -xvf BenElechi1.tar.gz
tar -xvf t3dh.tar.gz
tar -xvf F2.tar.gz
tar -xvf consph.tar.gz
tar -xvf SiO2.tar.gz
tar -xvf torso1.tar.gz
tar -xvf dielFilterV3real.tar.gz
tar -xvf RM07R.tar.gz
tar -xvf m_t1.tar.gz
tar -xvf crankseg_2.tar.gz
tar -xvf nd24k.tar.gz
tar -xvf TSOPF_RS_b2383.tar.gz
tar -xvf mouse_gene.tar.gz
tar -xvf human_gene1.tar.gz

echo "Expanding symmetric matrices"

../utilities/Matrix_Expander -i /parabolic_fem/parabolic_fem.mtx -o parabolic_fem_expanded.mtx
../utilities/Matrix_Expander -i /CurlCurl_3/CurlCurl_3.mtx -o CurlCurl_3_expanded.mtx
../utilities/Matrix_Expander -i /offshore/offshore.mtx -o offshore_expanded.mtx
../utilities/Matrix_Expander -i /nlpkkt80/nlpkkt80.mtx -o nlpkkt80_expanded.mtx
../utilities/Matrix_Expander -i /CO/CO.mtx -o CO_expanded.mtx
../utilities/Matrix_Expander -i /gsm_106857/gsm_106857.mtx -o gsm_106857_expanded.mtx
../utilities/Matrix_Expander -i /msdoor/msdoor.mtx -o msdoor_expanded.mtx
../utilities/Matrix_Expander -i /bmw3_2/bmw3_2.mtx -o bmw3_2_expanded.mtx
../utilities/Matrix_Expander -i /BenElechi1/BenElechi1.mtx -o BenElechi1_expanded.mtx
../utilities/Matrix_Expander -i /t3dh/t3dh.mtx -o t3dh_expanded.mtx
../utilities/Matrix_Expander -i /F2/F2.mtx -o F2_expanded.mtx
../utilities/Matrix_Expander -i /consph/consph.mtx -o consph_expanded.mtx
../utilities/Matrix_Expander -i /SiO2/SiO2.mtx -o SiO2_expanded.mtx
../utilities/Matrix_Expander -i /torso1/torso1.mtx -o torso1_expanded.mtx
../utilities/Matrix_Expander -i /dielFilterV3real/dielFilterV3real.mtx -o dielFilterV3real_expanded.mtx
../utilities/Matrix_Expander -i /m_t1/m_t1.mtx -o m_t1_expanded.mtx
../utilities/Matrix_Expander -i /crankseg_2/crankseg_2.mtx -o crankseg_2_expanded.mtx
../utilities/Matrix_Expander -i /nd24k/nd24k.mtx -o nd24k_expanded.mtx
../utilities/Matrix_Expander -i /mouse_gene/mouse_gene.mtx -o mouse_gene_expanded.mtx
../utilities/Matrix_Expander -i /human_gene1/human_gene1.mtx -o human_gene1_expanded.mtx

echo "Done preparing benchmark matrices."
