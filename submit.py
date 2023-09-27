import os, sys
from optparse import OptionParser

vs = {}
basenames = {}

usage = "usage: %prog [options]"
parser = OptionParser(usage, add_help_option=True)

parser.add_option("-a", "--nodeage", type="str", dest="nodeage", default=160,
                    help="New or old node template to edit")
parser.add_option("-c", "--numcpus", type="str", dest="numcpus", default=1,
                    help="Number of cpus/node working on the data set.")
parser.add_option("-n", "--number", type="str", dest="nM", default=40,
                    help="Number of points in each chunk (millions)")
parser.add_option("-m", "--mjd", type="str", dest="mjd", default=".",
                    help="MJD of data to search")
parser.add_option("-z", "--zmax", type="str", dest="zmax", default=160,
                    help="Maximum fourier drift (bins) to search")

(options, args) = parser.parse_args()

template = open("quick_template_%s.txt" % options.nodeage).read()

def submit(options):
    mjd, zmax, numchunks, numcpus = options.mjd, options.zmax, options.nM, options.numcpus
    print(f"Processing: {mjd}")
    job = f"Ter5All_new_files_{mjd}"
    out = template.replace("XXXX", mjd).replace("YYYY", zmax).replace("WWWW", numcpus).replace("ZZZZ", numchunks)
    outfile = open(job+".sh", "w")
    outfile.write(out)
    outfile.close()
    os.system("qsub %s.sh" % job)

# Submit
submit(options)