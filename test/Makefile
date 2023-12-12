all:
	$(RM) -rvf testwd;
	mkdir -pv testwd/logs;
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 1]: EM, then AMOVA with 1 levels in metadata, print distance matrix\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf --printDistanceMatrix 1 --minInd 2 -doDist 1 --maxEmIter 100 --em-tole 1e-10 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1 2> testwd/logs/test_s9_d1_1K_mtd1.log;
	bash -c "diff testwd/test_s9_d1_1K_mtd1.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 2]: Read the data from distance matrix file created in the previous step, then AMOVA with 1 levels in metadata\n";
	../ngsAMOVA -doEM 0 -doAMOVA 1 --in-dm testwd/test_s9_d1_1K_mtd1.distance_matrix.csv --minInd 2 -doDist 3 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_indm 2> testwd/logs/test_s9_d1_1K_mtd1_indm.log;
	bash -c "diff testwd/test_s9_d1_1K_mtd1_indm.amova.csv reference/test_s9_d1_1K_mtd1_indm.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 3]: EM, then AMOVA with 2 levels in metadata\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf --minInd 2 -doDist 1 --maxEmIter 100 --em-tole 1e-10 -m metadata_with_header_2lvl.tsv -f "Individual~Region/Population" -out testwd/test_s9_d1_1K_mtd2 2> testwd/logs/test_s9_d1_1K_mtd2.log;
	bash -c "diff testwd/test_s9_d1_1K_mtd2.amova.csv reference/test_s9_d1_1K_mtd2.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 4]: Read the data from distance matrix file, then do 1-level AMOVA using formula from metadata with 2 levels\n";
	../ngsAMOVA -doEM 0 -doAMOVA 1 --in-dm test_s9_d1_1K_mtd1_maxIter100_tole10.distance_matrix.csv --minInd 2 -doDist 3 --nThreads 0 -m metadata_with_header_2lvl.tsv -f "Individual~Population" -out testwd/indm_test_s9_d1_1K_mtd1 2> testwd/logs/indm_test_s9_d1_1K_mtd1.log;
	bash -c "diff testwd/indm_test_s9_d1_1K_mtd1.amova.csv testwd/test_s9_d1_1K_mtd1.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 5]: Use genotypes in VCF file to perform AMOVA\n";
	../ngsAMOVA -doEM 0 -doAMOVA 1 --in-vcf test_s9_d1_1K.vcf --minInd 2 -doDist 2 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_gt 2> testwd/logs/test_s9_d1_1K_mtd1_gt.log;
	bash -c "diff testwd/test_s9_d1_1K_mtd1_gt.amova.csv reference/test_s9_d1_1K_mtd1_gt.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 6]: Same as test 5, but use bcf input with more data and specify regions to use to get the same results as test 5\n";
	../ngsAMOVA -doEM 0 -doAMOVA 1 --in-vcf test_s9_d1_1K_extra.bcf --minInd 2 -doDist 2 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_mtd1_gt --region "chr22" 2> testwd/logs/test_s9_d1_1K_mtd1_gt.log;
	bash -c "diff testwd/test_s9_d1_1K_mtd1_gt.amova.csv reference/test_s9_d1_1K_mtd1_gt.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "[TEST 7]: 40 individuals, multithreaded EM\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 -doDist 1 --maxEmIter 50 --em-tole 1e-5 -P 4 --in-vcf sim_demes_v2-model1-1-rep0-d2.bcf -m sim_demes_v2-model1_metadata_2lvl_with_header.tsv -f "Individual~Region/Population" -o testwd/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5 2> testwd/logs/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.log;
	bash -c "diff testwd/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.amova.csv reference/sim_demes_v2-model1-1-rep0-d2_maxIter50tole5.amova.csv";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 8]: Neighbor Joining with individuals as tips. Use vcf individual names as tip labels\n";
	../ngsAMOVA -doDist 2 --in-vcf test_s9_d1_1K.vcf -f "Individual~Region/Population" -doPhylo 1 -o testwd/test_s9_d1_1K_mtd2_nj 2> testwd/logs/test_s9_d1_1K_mtd2_nj.log;
	bash -c "diff testwd/test_s9_d1_1K_mtd2_nj.newick reference/test_s9_d1_1K_mtd2_nj.newick";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 9]: Same as TEST1 but with same sites divided into 2 contigs\n";
	../ngsAMOVA -doEM 1 -doAMOVA 1 --in-vcf test_s9_d1_1K_2contigs.vcf --printDistanceMatrix 1 --minInd 2 -doDist 1 --maxEmIter 100 --em-tole 1e-10 -m metadata_with_header_1lvl.tsv -f "Individual~Population" -out testwd/test_s9_d1_1K_2contigs_mtd1 2> testwd/logs/test_s9_d1_1K_2contigs_mtd1.log;
	bash -c "diff testwd/test_s9_d1_1K_2contigs_mtd1.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
	@printf "\n\n\n";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 10]: Test if reading 3 gls and 10 gls works the same with the same data\n";
	../ngsAMOVA --in-vcf test_s9_d1_1K_3gls.vcf -doamova 1 -dodist 1 -doem 1 --formula "Individual~Population" -m metadata_with_header_1lvl.tsv --maxemiter 100 --em-tole 1e-10 -o testwd/test_s9_d1_1K_3gls_mtd1 2> testwd/logs/test_s9_d1_1K_3gls_mtd1.log;
	bash -c "diff testwd/test_s9_d1_1K_3gls_mtd1.amova.csv reference/test_s9_d1_1K_mtd1.amova.csv";
	@printf "\n\n\n";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 1]: Do ibd with ibdseq method implementation using default ibdseq params to ensure we get the same results as ibdseq\n";
	../ngsAMOVA --in-vcf test1_ibd_100.vcf --issim 1 -doibd 2 > testwd/test1_ibd_100_doIbd2.ibd;
	bash -c "diff testwd/test1_ibd_100_doIbd2.ibd reference/test1_ibd_100_doIbd2.ibd";
#
	@printf "\n\n\n";
	@printf "===============================================================================\n";
	@printf "\n\n\n[TEST 2]: Do ibd with both gl and ibdseq gt method implementations using input file where gls are simulated with an error rate of 0. Compare the gl and gt method (with errorprop0) to ensure we get the same results";
	../ngsAMOVA --in-vcf test2_ibd_1000_err0.vcf --issim 1 -doibd 2 --ibdseq-errorprop 0 > testwd/test2_ibd_1000_doIbd2_errorprop0.ibd
	../ngsAMOVA --in-vcf test2_ibd_1000_err0.vcf --issim 1 -doibd 1 > testwd/test2_ibd_1000_doIbd1.ibd
	bash -c "diff testwd/test2_ibd_1000_doIbd1.ibd testwd/test2_ibd_1000_doIbd2_errorprop0.ibd";
	bash -c "diff reference/test2_ibd_1000_doIbd1.ibd testwd/test2_ibd_1000_doIbd1.ibd";


