import os, sys
from optparse import OptionParser

basecommand = "pbsdsh -- python ${OUTDIR}/BBBB --numcpus ${NUMCPUS} --workdir ${OUTDIR} --zmax ${ZMAX} --flo ${FLO} --sigma ${SIGMA} --outdir ${OUTDIR} --numharm ${NUMHARM} --mjd ${MJD}"

def submit(options):
    mjd, zmax, numcpus, jobname, nodeage, searchmode = options.mjd, options.zmax, options.numcpus, options.jobname, options.nodeage, options.searchmode

    if nodeage == 'new':
        request = f"PBS -l nodes=5:new:ppn={numcpus}"
    elif nodeage == 'old':
        request = f"PBS -l nodes=7:old:ppn={numcpus}"
    else:
        print("Node age must be 'old' or 'new'!")
        sys.exit()
    
    print(f"Processing day: {mjd}")
    out = template.replace("XXXX", jobname).replace("WWWW", request).replace("AAAA", mjd).replace("EEEE", numcpus)
    if searchmode == 'chunk':
        searchscript = "/users/cchoza/chunk_search/short_analysis_Ter5All.py"
        searchcommand = basecommand.replace("BBBB", "short_analysis_Ter5All.py") + '  --number ${NUMCHUNK}'
        out = out.replace("ZZZZ", options.nM).replace("YYYY", zmax).replace("CCCC", searchscript).replace("DDDD", searchcommand)
    elif searchmode == 'full':
        searchscript = "/users/cchoza/chunk_search/full_analysis.py"
        searchcommand = basecommand.replace("BBBB", "full_analysis.py") + '  --fft'
        out = out.replace("YYYY", zmax).replace("CCCC", searchscript).replace("DDDD", searchcommand)
    else:
        print("Search mode must be 'full' or 'chunk'!")
        sys.exit()

    outfile = open(jobname+".sh", "w")
    outfile.write(out)
    outfile.close()
    os.system("qsub %s.sh" % jobname)


if __name__ == "__main__":
    vs = {}
    basenames = {}

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, add_help_option=True)

    parser.add_option("-a", "--nodeage", type="str", dest="nodeage", default=160,
                        help="New or old node template to edit")
    parser.add_option("-c", "--numcpus", type="str", dest="numcpus", default=1,
                        help="Number of cpus/node working on the data set.")
    parser.add_option("-n", "--number", type="str", dest="nM", default=40,
                        help="Number of points in each chunk (millions), for a chunk search.")
    parser.add_option("-m", "--mjd", type="str", dest="mjd", default=".",
                        help="MJD of data to search.")
    parser.add_option("-z", "--zmax", type="str", dest="zmax", default=160,
                        help="Maximum fourier drift (bins) to search, for acceleration searches.")
    parser.add_option("-j", "--jobname", type="str", dest="jobname", default=160,
                        help="Name for the job.")
    parser.add_option("-s", "--search", type="str", dest="searchmode", default=160,
                        help="Mode for the search. Can be 'short' for a short-chunk search, 'full' for a standard acceleration search.")

    (options, args) = parser.parse_args()
    template = open("search_template.txt").read()
    # Submit
    submit(options)