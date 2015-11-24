
#load testdata

myexp1 = read.table("exp1.txt")
myexp2 = read.table("exp2.txt")
mylist1=c('gene1','gene2','gene3','gene10')

#execute
mynvt <- NVTinit(mylist1,myexp1,myexp2,"N")

mynorm <- NVTnormalize(mynvt)

NVTplot(mynorm)
