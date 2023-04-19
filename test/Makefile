all:
	$(RM) -rvf testwd;
	mkdir -pv testwd;
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 1]: EM, then AMOVA with 1 levels in metadata, print distance matrix\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf -printMatrix 1 -isSim 1 -minInd 2 -doDist 1 -maxIter 100 -tole 1e-10 -m metadata_with_header_1lvl.tsv  -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1;
	bash -c "diff testwd/test_s9_d1_1K_mtd1.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 2]: Read the data from distance matrix file, then AMOVA with 1 levels in metadata\n";
	../ngsAMOVA -doEM 0 -doAMOVA 1 --in-dm testwd/test_s9_d1_1K_mtd1.distance_matrix.csv -isSim 1 -minInd 2 -doDist 1 -sqDist 0 -m metadata_with_header_1lvl.tsv  -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_indm;
	bash -c "diff testwd/test_s9_d1_1K_mtd1_indm.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 3]: EM, then AMOVA with 2 levels in metadata\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf -isSim 1 -minInd 2 -doDist 1 -maxIter 100 -tole 1e-10  -m metadata_with_header_2lvl.tsv  -f "Individual~Region/Population" -out testwd/test_s9_d1_1K_mtd2;
	bash -c "diff testwd/test_s9_d1_1K_mtd2.amova.csv reference/test_s9_d1_1K_mtd2.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 4]: Read the data from distance matrix file, then do 1-level AMOVA using formula from metadata with 2 levels\n";
	../ngsAMOVA -doEM 0 -doAMOVA 1 --in-dm test_s9_d1_1K_mtd1_maxIter100_tole10.distance_matrix.csv -isSim 1 -minInd 2 -doDist 1 -nThreads 0 -m metadata_with_header_2lvl.tsv  -sqDist 0  -f "Individual~Population" -out testwd/indm_test_s9_d1_1K_mtd1;
	bash -c "diff testwd/indm_test_s9_d1_1K_mtd1.amova.csv testwd/test_s9_d1_1K_mtd1.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 5]: Use genotypes in VCF file to perform AMOVA\n";
	../ngsAMOVA -doEM 0 -doAMOVA 2 --in-vcf test_s9_d1_1K.vcf -isSim 1 -minInd 2 -doDist 1 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_gt;
	bash -c "diff testwd/test_s9_d1_1K_mtd1_gt.amova.csv reference/test_s9_d1_1K_mtd1_gt.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 6]: Same as test 5, but use bcf input with more data and specify regions to use to get the same results as test 5\n";
	../ngsAMOVA -doEM 0 -doAMOVA 2 --in-vcf test_s9_d1_1K_extra.bcf -isSim 1 -minInd 2 -doDist 1 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_gt --region "chr22";
	bash -c "diff testwd/test_s9_d1_1K_mtd1_gt.amova.csv reference/test_s9_d1_1K_mtd1_gt.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "[TEST 7]: 40 individuals, multithreaded EM\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 -doDist 1 -maxIter 50 -tole 1e-5 -P 4 --in-vcf sim_demes_v2-model1-1-rep0-d2.bcf -m sim_demes_v2-model1_metadata_2lvl_with_header.tsv -f "Individual~Region/Population" -o testwd/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5;
	bash -c "diff testwd/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.amova.csv reference/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 8]: Neighbor Joining with individuals as tips\n";
	../ngsAMOVA -doAMOVA 2 -doDist 1 --in-vcf test_s9_d1_1K.vcf -m metadata_with_header_2lvl.tsv -f "Individual~Region/Population" -doNJ 1 -o testwd/test_s9_d1_1K_mtd2_nj;
	bash -c "diff testwd/test_s9_d1_1K_mtd2_nj.newick reference/test_s9_d1_1K_mtd2_nj.newick";
	@printf "\n\n\n";
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 9]: Same as TEST1 but with ancder file"
	../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf -printMatrix 1 -isSim 1 -minInd 2 -doDist 1 -maxIter 100 -tole 1e-10 -m metadata_with_header_1lvl.tsv  -f "Individual~Population" --ancDerFile test_s9_d1_1K_ancder.tsv -out testwd/test_s9_d1_1K_mtd1_ancder;
	bash -c "diff testwd/test_s9_d1_1K_mtd1_ancder.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";