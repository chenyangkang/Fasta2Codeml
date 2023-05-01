
          seqfile = Simple_multi_test.step7.phy            * Path to the alignment file
         treefile = Simple_multi_test.step7.tre           * Path to the tree file
          outfile = Simple_multi_test.null.res.txt            * Path to the output file

            noisy = 3              * How much rubbish on the screen
          verbose = 1              * More or less detailed report

          seqtype = 1              * Data type
            ndata = 1           * Number of data sets or loci
            icode = 0              * Genetic code 
        cleandata = 0              * Remove sites with ambiguity data?

            model = 2         * Models for ω varying across lineages
          NSsites = 2          * Models for ω varying across sites
        CodonFreq = 7        * Codon frequencies
          estFreq = 0        * Use observed freqs or estimate freqs by ML
            clock = 0          * Clock model
        fix_omega = 1         * Estimate or fix omega
            omega = 1        * Initial or fixed omega

    