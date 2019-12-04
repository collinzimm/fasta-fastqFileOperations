import sys
import matplotlib.pyplot as plt
import numpy as np

def readFile(file):
    file_handle = open(file, 'r')
    startsWithVar = '>'
    seqName = 'Default'
    desc = 'Default'
    sequence = ""
    seqList = []
    seqNameList = []
    seqDescList = []
    seq = file_handle.readlines()
    if file_handle.name.endswith('q'):
        temp = seq
        seq = []
        for i in temp:
            i = i.replace('@', '>')
            seq.append(i)
    for i in seq: #get name and desc
        if i.startswith('>') and len(i) < 80:
            seqName = i.split(' ')[0][1:]
            start = i.find(' ')
            desc = i[start:-1]
            seqList.append(sequence)
            seqNameList.append(seqName)
            seqDescList.append(desc)
            sequence = ""
        else: #seq
            if i.startswith('+') == False and i.count('/') == 0 and i.count('>') == 0:
                sequence += i.strip()
    seqList.append(sequence)
    if seqList.count(''):
        seqList.remove('')
    return seqName, desc, seqList;


#prints in fasta format
def printInFasta(SeqName, lines, SeqDesc, n = 100):
    seqBase = ''
    for i in lines:
        if i.startswith('>'):
            seqName = i.split(' ')[0][1:]
        else:
            seqBase = seqBase+i.strip()

    for i in range(0, len(seqBase), n):
        line = seqBase[i:(i+n)]
        for j in range(0, len(line), 10):
            tenSeq = line[j:(j+10)]
            print(tenSeq, end = '')
        print()


#same as printinFasta but with headers, seems redundent
def printWithRuler(lines, delim, n=100):
    # print headers
    if lineNum:
        print('     ', end = '')
        for i in range(int(n/10)):
            print('         ' + str(i + 1), end = delim)
        print('\n')
        print('Line ', end = '')
        for j in range(int(n/10)):
            print('1234567890', end = delim)
        print('')
    seqBase = ''
    #split the lines
    for i in lines:
        if i.startswith('>'):
            seqName = i.split(' ')[0][1:]
        else:
            seqBase = seqBase+i.strip()
    #print the lines
    for i in range(0, len(seqBase), n):
        line = seqBase[i:(i+n)]
        num = '   '+str(int(i/n)+1)+' '
        print(num, end = '')
        for j in range(0, len(line), 10):
            tenSeq = line[j:(j+10)]
            print(tenSeq, end = delim)
        print()


#return gc count percent as decimal
def gcContent(lines):
    totalLength = 0
    for i in lines:
        if i.count('>') == 0:
            totalLength = totalLength + len(i)
            countC = i.count('C')
            countG = i.count('G')
    result = (countC + countG)/totalLength
    return result


#same as gcCount but with more chars,
def nucleotideCounter(lines):
    countA = 0
    countT = 0
    countG = 0
    countC = 0
    countN = 0
    totalLength = 0
    for i in lines:
        if i.count('>') == 0:
            totalLength = totalLength + len(i)
            countA = countA + i.count('C')
            countT = countT + i.count('T')
            countG = countG + i.count('G')
            countC = countC + i.count('C')
            countN = countN + i.count('N')
    countO = totalLength - len(lines)
    return countA, countT, countG, countC, countN, countO;

#DNA to Protein
def translation(seq):
    table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
    protein = ''
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon.count('N') == 0:
                protein = protein + table[codon]
    return protein


def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    revComp = "".join(complement.get(base, base) for base in reversed(seq))
    return revComp

def ORF(temp):
    temp = temp.replace('N', '') #dna seq
    orfList1 = []
    orfList2 = []
    orfList3 = []
    orfList4 = []
    orfList5 = []
    orfList6 = []
    orfListList = []
    mTrack = False
    endTrack = False
    frames = [1,2,3,-1,-2,-3]
    for num in frames:
        protein =[]
        dna = temp #reassign dna to origional seq
        #prep dna seq frame for translation
        if num == 1:
            dna = dna
        if num == 2:
            dna = dna[1:]
        if num == 3:
            dna = dna[2:]
        #reverse compliment
        if num == -1:
            dna = reverseComplement(dna)
        if num == -2:
            dna = reverseComplement(dna)
            dna = dna[1:]
        if num == -3:
            dna = reverseComplement(dna)
            dna = dna[2:]
        #translation
        for j in range(0, len(dna), 3):
                threeSeq = dna[j:(j+3)]
                protein += translation(threeSeq)
        orf = []
        #find orf
        for i in range(0, len(protein), 1):
            if protein[i] == "M": #find start
                mTrack = True
            if protein[i] == "*":#ends the orf
                endTrack = True
                orf.append('*')
            if mTrack is True and endTrack is False: #adds to orf when M is found
                orf.append(protein[i])
            if endTrack is True: #validates orf when end is found
                if len(orf) > 10: #checks length of orf
                    if num == 1:
                        orfList1.append(orf) #adds orf to list
                    if num == 2:
                        orfList2.append(orf) #adds orf to list
                    if num == 3:
                        orfList3.append(orf) #adds orf to list
                    if num == -1:
                        orfList4.append(orf) #adds orf to list
                    if num == -2:
                        orfList5.append(orf) #adds orf to list
                    if num == -3:
                        orfList6.append(orf) #adds orf to list
                orf = [] #clears orf
                endTrack = False
                mTrack = False
                orfList1.sort(key = len, reverse= True)
                orfList2.sort(key = len, reverse= True)
                orfList3.sort(key = len, reverse= True)
                orfList4.sort(key = len, reverse= True)
                orfList5.sort(key = len, reverse= True)
                orfList6.sort(key = len, reverse= True)
        orfListList.append(orfList1)
        orfListList.append(orfList2)
        orfListList.append(orfList3)
        orfListList.append(orfList4)
        orfListList.append(orfList5)
        orfListList.append(orfList6)
    lists = [1,2,3,4,5,6]
    orfNumber = ["+1", "+2", "+3", "-1", "-2", "-3"]
    for num in lists:
        orfList = orfListList[num-1]
        if len(orfList) > 0:
            for k in range(len(orfList)):
                print("ORF ("+ str(orfNumber[num-1]) +"): ", end ="")
                for x in range(len(orfList[k])):
                    print(orfList[k][x], end = "")
                print()
        else:
            print("ORF ("+ str(orfNumber[num-1]) +"): ", end ="")
            print("NO VALID ORF")


def gcPlotCse(dnaList, binNum):
    for i in range(len(dnaList)):
        tenthList = []
        xAxis = []
        yAxis = []
        temp = dnaList[i]
        dna = temp
        tenth = int(len(dna)/10)
        for k in range(1,11,1):
            tenthList.append(tenth * k)
        dnaTenths = []
        begin = 0
        for k in range(len(tenthList)):
            dnaTenths.append(dna[begin:tenthList[k]])
            begin = tenthList[k]
            tempList = [dnaTenths[k]]
            xAxis.append(gcContent(tempList))
    if len(dnaList) == 1:
        fig, (ax1) = plt.subplots(1, 1) #row, col
        ax1.hist(xAxis, binNum)
        ax1.title.set_text("Seq 1")
        ax1.set_xlabel("GC content")
        ax1.set_ylabel("Counts")
        plt.show()
    if len(dnaList) == 2:
        fig, (ax1, ax2) = plt.subplots(1, 2) #row, col
        ax1.hist(xAxis, binNum)
        ax1.title.set_text("Seq 1")
        ax1.set_xlabel("GC content")
        ax1.set_ylabel("Counts")
        ax2.hist(xAxis, binNum)
        ax2.title.set_text("Seq 2")
        ax2.set_xlabel("GC content")
        ax2.set_ylabel("Counts")
        plt.show()
    if len(dnaList) == 3:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3) #row, col
        ax1.hist(xAxis, binNum)
        ax1.title.set_text("Seq 1")
        ax1.set_xlabel("GC content")
        ax1.set_ylabel("Counts")
        ax2.hist(xAxis, binNum)
        ax2.title.set_text("Seq 2")
        ax2.set_xlabel("GC content")
        ax2.set_ylabel("Counts")
        ax3.hist(xAxis, binNum)
        ax3.title.set_text("Seq 3")
        ax3.set_xlabel("GC content")
        ax3.set_ylabel("Counts")
        plt.show()
    if len(dnaList) >= 4:
        fig, axs = plt.subplots(2,2) #row, col
        axs[0, 0].hist(xAxis, binNum)
        axs[0, 0].title.set_text("Seq 1")
        axs[0, 0].set_xlabel("GC content")
        axs[0, 0].set_ylabel("Counts")
        axs[0, 1].hist(xAxis, binNum)
        axs[0, 1].title.set_text("Seq 2")
        axs[0, 1].set_xlabel("GC content")
        axs[0, 1].set_ylabel("Counts")
        axs[1, 0].hist(xAxis, binNum)
        axs[1, 0].title.set_text("Seq 3")
        axs[1, 0].set_xlabel("GC content")
        axs[1, 0].set_ylabel("Counts")
        axs[1, 1].hist(xAxis, binNum)
        axs[1, 1].title.set_text("Seq 4")
        axs[1, 1].set_xlabel("GC content")
        axs[1, 1].set_ylabel("Counts")
        plt.show()

#end of functions

file = sys.argv[1]
seqNameList, seqDescList, linesList = readFile(file)
seqExamNum = 1 #test number to examine given by user
print("There are " + str(len(linesList)) + " sequence(s) in the file: " + file)
seqNumList = list(range(1, len(linesList)+1))
seqNumStr = "["
for i in seqNumList:
    if i != len(seqNumList):
        seqNumStr = seqNumStr + str(seqNumList[i-1]) + "|"
    else:
        seqNumStr = seqNumStr + str(seqNumList[i-1]) + "]"
n = 100
seqBase =''
Q1 = ''
Q2 = ''
lineNum = False
validAns= False
# format selection
while Q1 != 'N' and Q1 != 'Y' and validAns == False:

    while validAns == False:
        n = input('How many nucleotides per line [50|100]?')
        if n == '50' or n == '100':
            validAns = True
            n = int(n)
    validAns = False
    while validAns == False:
        seqExamNum = input('Which sequence do you want to examine ' + seqNumStr +'?')
        for i in seqNumList:
            if seqExamNum == str(seqNumList[i-1]):
                validAns = True

    Q1 = input('FASTA format? (Y/N)')
    if Q1 == 'Y':
        lineNum = False
    if Q1 == 'N':
        lineNum = True

    # space or no space
    if Q1 == 'N':
        while Q2 != 'N' and Q2 != 'Y':
            Q2 = input('Space delimiter? (Y/N)')
            if Q2 == 'Y':
                s = ' '
            if Q2 == 'N':
                s = ''

seqExam = int(seqExamNum) - 1
lines = linesList[seqExam]
seqName = seqNameList[seqExam]
seqDesc = seqDescList[seqExam]

#call print function
if lineNum:
    printWithRuler(lines, s, n)
else:
    printInFasta(seqName, lines, seqDesc, n)
var = -8
countA, countT, countG, countC, countN, countO = nucleotideCounter(lines)
print('A = ' + str(countA))
print('T = ' + str(countT))
print('G = ' + str(countG))
print('C = ' + str(countC))
print('N = ' + str(countN))
print('Other = ' + str(countO))
tempLine = [lines]
print('GC content = ' + str(gcContent(tempLine)))

#prep list for gcPlot here

for i in lines:
    if i.startswith('>'):
        seqName = i.split(' ')[0][1:]
    else:
        seqBase = seqBase+i.strip()

ORF(seqBase)
n = 5
gcPlotCse(linesList, n)


