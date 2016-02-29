import fileinput

for i,line in enumerate(fileinput.input()):
    print line[-200:-1]
