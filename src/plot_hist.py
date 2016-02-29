import matplotlib.pyplot as plt
import cPickle as pickle
import fileinput
import os
import math

def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.

    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.

    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.

    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.

    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.

    """
    
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
 
    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
 
    # The final path to save to
    savepath = os.path.join(directory, filename)
 
    if verbose:
        print("Saving figure to '%s'..." % savepath),
 
    # Actually save the figure
    plt.savefig(savepath)
    
    # Close it
    if close:
        plt.close()
 
    if verbose:
        print("Done")
######################################
        
def main():
    rna_out_total = []   
    for line in fileinput.input():
        line = line.rstrip()
        if float(line)>0:
            rna_out_total.append(math.log(float(line)+0.001))
             
    plt.hist(rna_out_total)
    plt.title("RNASeq Median Expression of 12 MCL Patients")
    plt.xlabel("ln(RPKM)")
    plt.ylabel("Frequency")
    save('MCLRNASeq_out_all','png')
main()