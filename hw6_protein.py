"""
Protein Sequencing Project
Name:
Roll Number:
"""

from tkinter import Label
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file = open(filename, "r").read()
    split = file.splitlines()
    string = "".join(split)
    return string


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    rnalist = []
    stoplist = ["UAA", "UGA", "UAG"]
    dna1 = dna[startIndex:]
    dna2 = dna1.replace("T", "U")

    for i in range(0,len(dna2),3):
        rn = dna2[i:i+3]
        if rn in stoplist:
            rnalist.append(rn)
            break
        rnalist.append(rn)

    return rnalist


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    dic = {}
    import json
    f = open(filename)
    data = json.load(f)
    for key,values in data.items():
        for i in values:
            re= i.replace("T", "U")
            dic[re] = key
    return dic


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protlist = []
    for i in codons:
        if i == "AUG" and len(protlist) == 0:
            protlist.append("Start")
        else:
            protlist.append(codonD[i])
    return protlist


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    synlist = []
    x = readFile(dnaFilename)
    y = makeCodonDictionary(codonFilename)
    i = 0
    unused = 0

    while i < len(x):
        sli = x[i:i+3]
        if sli == "ATG":
            rna = dnaToRna(x, i)
            prolist = generateProtein(rna, y)
            synlist.append(prolist)
            i += 3 * len(prolist)
        else:
            i += 1
            unused += 1
    print("total number of bases" , len(x), "unused base count", unused, "total number of protiens", len(synlist) )
    # print(synlist)
    return synlist


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    finallist = []
    for i in proteinList1:
        if i in proteinList2:
            finallist.append(i)
    return finallist


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    final = []
    for i in proteinList:
        for j in i:
            final.append(j)
    return final


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    final = {}
    for i in aaList:
        count = aaList.count(i)
        final[i] = count
    return final


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    list1 = combineProteins(proteinList1)
    list2 = combineProteins(proteinList2)
    dict1 = aminoAcidDictionary(list1)
    dict2 = aminoAcidDictionary(list2)
    
    l = []
    
    for i in dict1:
        if i not in l and i != "Start" and i != "Stop":
            l.append(i)
    for i in dict2:
        if i not in l and i != "Start" and i != "Stop":
            l.append(i)    
   
    l2 = []     
    for amo in l:  
        f1 = 0
        f2 = 0
        if amo in list1:
            f1 = (dict1[amo]/ len(list1))
        if amo in list2:
            f2 = (dict2[amo]/len(list2))
        dif = f2 -f1
        if (dif > cutoff) or (dif < -cutoff):
            l2.append([amo, f1, f2])
    return l2


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    l1 = []
    for i in commonalities:
        i.remove("Start")
        i.remove("Stop")    
        if len(i) > 1:
            j = "-".join(i)
            l1.append([j])
        else:
            if i not in l1:
                l1.append(i)
    l2 = sorted(l1)
    print("The following proteins occurred in both DNA Sequences:")
    for a in l2:
        for b in a:
            print(b)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for val in differences:
        word = val[0]
        seq1 = round(val[1] *100,2)
        seq2 = round(val[2]*100,2)
        print(f"{word}: {seq1} % in seq1, {seq2} % in seq2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    list = []
    l1 = combineProteins(proteinList1)
    l2 = combineProteins(proteinList2)
    for i in l1:
        if i not in list:
            list.append(i)
    for j in l2:
        if j not in list:
            list.append(j)
    final = sorted(list)
    # print(final)
    return final


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    l1 = combineProteins(proteinList)
    dic1 = aminoAcidDictionary(l1)
    l2 = []
    for i in labels:
        if i in dic1:
            freq = dic1[i]/len(l1)
            l2.append(freq)
        else:
            l2.append(0)
    return l2


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList="black"):
    import matplotlib.pyplot as plt

    w = 0.35  # the width of the bars

    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1, edgecolor= edgeList)
    plt.bar(xLabels, freqList2, width= w, align='edge', label=label2, edgecolor= edgeList)

    plt.xticks(rotation="vertical")
    plt.legend()
    plt.title("Compare Human and Elephant")

    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    diff = []
    for i in biggestDiffs:
        diff.append(i[0])
    l1 = []
    for i in labels:
        if i in diff:
            l1.append("black")
        else:
            l1.append("white")
    return l1


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():

    proteins1 = synthesizeProteins("data/Human_p53.txt", "data/codon_table.json")
    proteins2 = synthesizeProteins("data/Elephant_p53.txt", "data/codon_table.json")
    common = commonProteins(proteins1, proteins2)
    amodif = findAminoAcidDifferences(proteins1,proteins2,0.005)
    disply = displayTextResults(common, amodif)
    labels = makeAminoAcidLabels(proteins1,proteins2)
    f1 = setupChartData(labels, proteins1)
    f2 = setupChartData(labels, proteins2)
    edges = makeEdgeList(labels, amodif)
    createChart(labels, f1, "Human", f2, "Elephant", edgeList=edges)
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()

    ## Uncomment these for Week 2 ##
    # test.testCommonProteins()
    # test.testCommonProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()
    # test.testCreateChart()
    # test.testMakeEdgeList()
    
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    # test.testMakeAminoAcidLabels()
    # test.testSetupChartData()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
