from utilities import *
import random
sample = list("dskfslbndlkbndkfbndlfbkdnflbkfbndlfkbndflb")

b = 'sfkjsd0jkl'
for i in range(1,1000000):
    a = "".join(random.sample(sample,10)) #'sfdk0fjsdlk'
    c = similar_comp(a,b)
print c