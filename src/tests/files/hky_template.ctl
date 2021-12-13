     seqfile = #SEQFILE
    treefile = #TREEFILE
     outfile = #OUTPUTFILE

       noisy = 1
     verbose = 0
     runmode = 0

       model = 4   * HKY model
       Mgene = 0   * No variation of parameters across genes
       ndata = 1   * How many datasets to analyze
       nhomo = 5   * Specify the number of model parameters and where those models are
   fix_kappa = 0   * One set of parameters estimated per branch
                   * 1 or 2 do different, weird things, see documentation

       clock = 0   * Rates are free to vary across branches (use an unrooted tree)
   fix_alpha = 1   * Fix alpha parameter 
       alpha = 0.  * Infinity (single rate for all sites) 
    
       getSE = 0   * Don't want estimates of standard errors
RateAncestor = 0   * Don't want that extra analysis
   cleandata = 1   * Sites involving ambiguity or alignment gaps are removed
      method = 0   * Can't use nhomo option and method
                   * So all parameters and branch lengths simultaneously
 fix_blength = 0   * Ignore branch lengths in tree file

