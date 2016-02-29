from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
figure = plt.figure()
venn2(subsets = (33, 200, 1))
plt.show()