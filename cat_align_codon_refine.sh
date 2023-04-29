for name in $(<unique_gene_names.txt)
do
{
	cat unique_cds_names.txt|grep ^$name. > ${name}.cds_list.txt
    
	echo '''

	cd /beegfs/store4/chenyangkang/06.ebird_data/40.Shift_and_Genome/02.genomes/genome_alignment

	for name in $(<'"${name}"'.cds_list.txt)
	do
	{
		cat cds_sequence_no_concat/*/${name}.fasta > all_multispecies_fasta/${name}.fasta
		wait
        
		filterout_full_N_seq.py all_multispecies_fasta/${name}.fasta
		wait
        
		muscle -in all_multispecies_fasta/${name}.fasta -out all_multispecies_aln/${name}.aln.fasta -maxiters 5 -diags
		wait
        
		java -jar -Xmx6G /beegfs/store4/chenyangkang/software/macse_v2.07.jar -prog refineAlignment \
		-align all_multispecies_aln/${name}.aln.fasta \
		-out_AA codon_AA_alignment/${name}.aln.AA.fasta \
		-out_NT codon_NT_alignment/${name}.aln.NT.fasta
		wait
        
		java -jar -Xmx6G \
		/beegfs/store4/chenyangkang/software/macse_v2.07.jar \
		-prog exportAlignment \
		-codonForInternalFS NNN \
		-codonForInternalStop NNN \
		-charForRemainingFS N \
		-codonForExternalFS NNN \
		-codonForFinalStop NNN \
		-align codon_NT_alignment/${name}.aln.NT.fasta \
		-out_NT codon_NT_alignment/${name}.aln.NT.Trimed.fasta
        
		#wait
		#raxml -f a -x 42 -p 42 -# 100 -m GTRGAMMA -s test.phy -n test -T 1
		## cannot run raxml here. raxml need to be run on gene bases
		wait
	}
	done
	wait
	
	rm '"${name}"'.cds_list.txt
	wait
    	sleep 1

	concatenate_alignment.py '"${name}"'
	wait
	sleep 1

	co_filter_aln_and_tree.py '"${name}"'
	wait
	sleep 1

	make_paml_ctl_file.py '"${name}"'
	wait
	sleep 1

	cd /beegfs/store4/chenyangkang/06.ebird_data/40.Shift_and_Genome/02.genomes/genome_alignment/run_PAML/'"${name}"'
	wait

	codeml '"${name}"'.alternative.ctl > alt.log
	wait
	sleep 1

	codeml '"${name}"'.null.ctl > null.log
	wait
	sleep 1

	parse_paml_output.py \
	--null_file '"${name}"'.null.res.txt \
	--alt_file '"${name}"'.alternative.res.txt \
	--out_file ./Stats_rest.csv \
	--df 1

	''' > ${name}.aln.sh
	chmod 755 ${name}.aln.sh

	while true
	do
	{
		myjobcount=`qstat |grep "chenyangkang"|wc -l`
		if [[ $myjobcount -lt 150 ]];then
			qsub -l nodes=1:ppn=1,mem=6G ${name}.aln.sh
			break
		fi
		sleep 5
	}
	done

	rm ${name}.aln.sh
	sleep 0.2

}
done
