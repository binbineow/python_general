#this file contains useful functions I wrote in different occasians
#one line description at the top
#import should be included inside of the function

#save matplotlib plot into a file
def save(path, ext='png', close=True, verbose=True):
    import os
    import matplotlib.pyplot as plt
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

#get a list from a file name###        
def get_list_from_file(file_name0):
    file0 = open(file_name0,'r')
    list_out = []
    for x in file0:
        x = x.rstrip()
        list_out.append(x)
    return list_out
######################################

def get_val_from_dict(dict0,make_set=True):
    list_val = []
    for _,val0 in dict0.iteritems():
        list_val0.append(val0)
    if make_set:
        list_val = list(set(list_val))
    return list_val
######################################