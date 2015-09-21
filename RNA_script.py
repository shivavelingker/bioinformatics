#-------------------------------------------------------------------------------
# Name:        Pre-mRNA Computational Script
# Purpose:     To retrieve exonic and intronic density information from GTF and Bedgraph files and compute pre-mRNA fractions as determined by various tissues
#
# Author:      Shiva Velingker, shivavelingker@student.lsmsa.edu
#              Louisiana School for Math, Science, and the Arts
# Created:     07/01/2014
# Copyright:   (c) Shiva Velingker 2014
#-------------------------------------------------------------------------------
import linecache, shutil, os, gzip, time

#RETRIEVE ALL EXONS FROM GTF FILE
def getKeyGTF(item):
    return [str(item[0]), int(item[4])]
def getExons(gtfFile, output, start):
    newOutput = str(output)+"\\"+str(start)+"new.txt"
    myFile = open(gtfFile, "r")
    newFile = open(newOutput, "w")
    geneFile = open(str(output)+"\\"+str(start)+"gene.txt", "w")
    array = []

    #skip any comments at the beginning of file
    line = myFile.readline().split('#')
    if(len(line)>1):
        while(len(line)>1):
            line = myFile.readline().split('#')
    line = line[0]

    while line: #until end of GTF
        line = line.split('\t')

        if(str(line[2])=="gene"): #if the line is detected to be exon
            geneID = line[8].split(';')[0]
            geneID = geneID.strip("gene_id ")
            geneID = geneID.strip("\"")

            start = line[3]
            stop = line[4]
            chrom = line[0]

            geneFile.write(str(chrom) +'\t'+ str(geneID) +'\t'+ str(int(start)) +'\t'+ str(int(stop)) +'\n')

        if(str(line[2])=="exon"): #if the line is detected to be exon
            geneID = line[8].split(';')[0]
            geneID = geneID.strip("gene_id ")
            geneID = geneID.strip("\"")
            geneName = line[8].split(';')[3]
            geneName = geneName.strip("gene_name ")
            geneName = geneName.strip("\"")
            biotype = line[8].split(';')[5]
            biotype = biotype.strip("gene_biotype ")
            biotype = biotype.strip("\"")

            chrom = line[0]
            start = line[3]
            stop = line[4]

            array.append([chrom,geneID,geneName,biotype,start,stop])

        line = myFile.readline()
    array = sorted(array, key=getKeyGTF)
    for counter in range(0, len(array)):
        newFile.write(str(array[counter][0]) +'\t'+ str(array[counter][1]) +'\t'+ str(array[counter][2]) +'\t'+ str(array[counter][3]) +'\t')
        newFile.write(str(int(array[counter][4])) +'\t'+ str(int(array[counter][5])) +'\n')
    myFile.close()
    newFile.close()
    geneFile.close()

#CONDENSE OVERLAP BETWEEN FILES
def overlap(first, second, lineNum, num_lines):
    first = first.split('\t')
    second = second.split('\t')

    if(lineNum == num_lines): #if end of file
        return [first[0], first[1], first[2], first[3], first[4], first[5], 'EXON']

    if(int(first[4])<=int(second[4]) and int(second[4])<=int(first[5])): #IS overlap
        if(str(first[1]) != str(second[1])):
            return [0,0,0,0,0, 'OWN']
        if int(first[4]) <= int(second[4]):
            newF = int(first[4])
        else:
            newF = int(second[4])
        if int(first[5]) <= int(second[5]):
            newE = int(second[5])
        else:
            newE = int(first[5])
        return [first[0], first[1], first[2], first[3], newF, newE, 'OWN']
    else: #no overlap
        return [first[0], first[1], first[2], first[3], first[4], first[5], 'EXON']

def condExons(output, start):
    loop = False
    num_lines = sum(1 for line in open(str(output)+"\\"+str(start)+"new.txt"))

    newFile = open(str(output)+"\\"+str(start)+"exons.txt", "w+")
    discardFile = open(str(output)+"\\"+str(start)+"discard.txt", "a")

    lineNum = 1
    while lineNum <= num_lines:
        thisLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum).strip('\n')
        nextLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum+1).strip('\n')

        writeLine = overlap(thisLine, nextLine, lineNum, num_lines) #finds is overlap exists between 2 exons

        if str(writeLine[0]) == '0': #2 exons come from different genes
            lineNum = lineNum + 2
            thisLine = thisLine.split('\t')
            nextLine = nextLine.split('\t')
            discardFile.write(str((thisLine[0])) +'\t'+ str(thisLine[1]) +'\t'+ str(thisLine[2])+'\t'+ str(thisLine[3]) +'\t'+ str(int(thisLine[4])) +'\t'+ str(int(thisLine[5])) +'\tDISCARDED EXON\n')
            discardFile.write(str((nextLine[0])) +'\t'+ str(nextLine[1]) +'\t'+ str(nextLine[2])+'\t'+ str(thisLine[3]) +'\t'+ str(int(nextLine[4])) +'\t'+ str(int(nextLine[5])) +'\tDISCARDED EXON\n')
        elif str(writeLine[6]) == 'OWN': #condense a pair of overlapping genes and moves to next pair
            lineNum = lineNum + 2
            loop = True #will tell program that there are still overlaps
            newFile.write(listToStr(writeLine) +'\n')
        else: #thisExon had no overlaps with the next line
            lineNum = lineNum + 1
            newFile.write(listToStr(writeLine) +'\n')

    newFile.close()
    discardFile.close()
    linecache.clearcache()

    shutil.copyfile(str(output)+"\\"+str(start)+"exons.txt", str(output)+"\\"+str(start)+"new.txt")
    return loop

#ADD INTRONS
def insert(first, second, lineNum, num_lines):
	first = first.split('\t')
	second = second.split('\t')

	if(lineNum == num_lines): #if end of line
		return [first[0], first[1], first[2], first[3], first[4], first[5], 'EXON']

	if(str(first[1]) != str(second[1])): #if different genes between 2 lines
		return [0,0,0,0,0,0,0]
	else: #return code to insert intron
		start = int(first[5])+1
		end = int(second[4])-1
		return [first[0], first[1], first[2], first[3], start, end, 'INTRON']

def getIntrons(output, start):
    num_lines = sum(1 for line in open(str(output)+"\\"+str(start)+"new.txt"))

    newFile = open(str(output)+"\\"+str(start)+"exons.txt", "w+")

    lineNum = 1
    while lineNum <= num_lines:
        thisLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum).strip('\n')
        nextLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum+1).strip('\n')

        writeLine = insert(thisLine, nextLine, lineNum, num_lines) #finds location of intron
        thisLine = thisLine.split('\t')
        if str(writeLine[6]) == '0': #difference in genes
            lineNum = lineNum + 1
            newFile.write(listToStr(thisLine) +'\n')
        elif str(writeLine[6]) == 'INTRON': #insert intron after exon
            lineNum = lineNum + 1
            newFile.write(listToStr(thisLine) +'\n')
            newFile.write(listToStr(writeLine) +'\n')
        else: #end of file
            lineNum = lineNum + 1
            newFile.write(listToStr(writeLine) +'\n')

    newFile.close()
    linecache.clearcache()

    shutil.copyfile(str(output)+"\\"+str(start)+"exons.txt", str(output)+"\\"+str(start)+"new.txt")

#ALIGN WITH BEDGRAPH
def over(bedLine, readLine):
    if(str(bedLine[0])==str(readLine[0]) and int(bedLine[1]) > int(readLine[5])): #if bedgraph doesn't support exon location
        return True
    elif(str(bedLine[0])>str(readLine[0])): #if bedgraph doesn't support chrom loc
        return True
    else:
        return False

def end(bedLine, readLine): #checks to see if exon/intron has reached end of bed alignment
    if(int(bedLine[2])>=int(readLine[5])):
        return True
    else:
        return False

def intersect(bedLine, readLine): #checks for alignment between bed and exon/intron
    if(str(bedLine[0])!=str(readLine[0])): #if not from same chromosomes
        return False
    if(int(bedLine[1])<=int(readLine[4]) and int(readLine[4])<=int(bedLine[2])): #b r b
        return True
    elif(int(readLine[4])<=int(bedLine[1]) and int(bedLine[1])<=int(readLine[5])): #r b r
        return True
    else: #no intersection
        return False

def densities(bedgraphFile, output, start):
    num_lines = sum(1 for line in open(str(output)+"\\"+str(start)+"exons.txt"))
    lineNum = 1
    
    readFile = open(str(output)+"\\"+str(start)+"exons.txt", "r")
    newFile = open(str(output)+"\\"+str(start)+"new.txt", "w")
    discardFile = open(str(output)+"\\"+str(start)+"discard.txt", "a")
    readLine = readFile.readline().strip('\n').split('\t')

    cont = True
    oldChrom = readLine[0]
    print "\tBeginning with chromosome: "+str(oldChrom)
    bedFile = gzip.open(bedgraphFile, "r")
    bedLine = bedFile.readline()

    exonTotal = 0
    intronTotal = 0
    exonDiscard = 0
    intronDiscard = 0

    bedLine = bedFile.readline().strip('\n').split('\t')
    while lineNum<=num_lines:
        exonDens = 0
        intronDens = 0
        begin = False
        while(not begin and lineNum<=num_lines): #first get to read-coverage point
            coverage = over(bedLine, readLine)
            if(coverage): #if exon isn't supported in bedgraph, skip it
                if(str(readLine[6])=='EXON'):
                    exonDiscard = exonDiscard + 1
                    discardFile.write(str(readLine[0]) +'\t'+ str(readLine[1]) +'\t'+ str(readLine[2]) +'\t'+ str(readLine[3]) +'\t'+ str(int(readLine[4])) +'\t'+ str(int(readLine[5])) +'\t'+ 'EXON NOT SUPPORTED BY BEDGRAPH\n')
                else:
                    intronDiscard = intronDiscard + 1
                    discardFile.write(str(readLine[0]) +'\t'+ str(readLine[1]) +'\t'+ str(readLine[2]) +'\t'+ str(readLine[3]) +'\t'+ str(int(readLine[4])) +'\t'+ str(int(readLine[5])) +'\t'+ 'INTRON NOT SUPPORTED BY BEDGRAPH\n')
                readLine = readFile.readline().strip('\n').split('\t')
                lineNum = lineNum + 1
            begin = intersect(bedLine, readLine) #find intersection with new/old exon
            if(not begin and not coverage): #if no intersection BUT exon is still supported, look to next bedLine
                bedLine = bedFile.readline()
                if(bedLine):
                    bedLine = bedLine.strip('\n').split('\t')
                else:
                    lineNum = num_lines+1

        if(lineNum<=num_lines):
            tooShort = end(bedLine, readLine)
            
            if(tooShort): #if E/I starts/stops within bedLine
                if(str(readLine[6]) == 'EXON'):
                    exonDens = exonDens + (int(readLine[5]) - int(readLine[4]) + 1) * int(bedLine[3])
                else:
                    intronDens = intronDens + (int(readLine[5]) - int(readLine[4]) + 1) * int(bedLine[3])
            else: #if not, get density for first alignment
                if(str(readLine[6]) == 'EXON'):
                    exonDens = exonDens + (int(bedLine[2]) - int(readLine[4]) + 1) * int(bedLine[3])
                else:
                    intronDens = intronDens + (int(bedLine[2]) - int(readLine[4]) + 1) * int(bedLine[3])

            if(not tooShort):
               bedLine = bedFile.readline().strip('\n').split('\t') #move to next bedgraph line
            else:
               begin = False

        while(begin and lineNum<=num_lines): #get densities for the rest
            begin = intersect(bedLine, readLine) #checks to see if exon/intron still aligns with bedgraph
            if(begin):
                if(end(bedLine, readLine)): #check to see if this bed's location terminates the exon/intron
                    if(readLine[6] == 'EXON'):
                        exonDens = exonDens + ((int(readLine[5]) - int(bedLine[1]) + 1)*int(bedLine[3]))
                    else:
                        intronDens = intronDens + ((int(readLine[5]) - int(bedLine[1]) + 1)*int(bedLine[3]))
                    begin = False
                else: #if no termination, calculate and move to next bedgraph line
                    if(readLine[6] == 'EXON'):
                        exonDens = exonDens + ((int(bedLine[2]) - int(bedLine[1]))*int(bedLine[3]))
                    else:
                        intronDens = intronDens + ((int(bedLine[2]) - int(bedLine[1]))*int(bedLine[3]))
                    bedLine = bedFile.readline().strip('\n').split('\t')
        if(lineNum<=num_lines and float(int(readLine[5])-int(readLine[4]) + 1) != 0):
            denom = float(int(readLine[5])-int(readLine[4]) + 1)
            eiSize = int(readLine[5])-int(readLine[4]) + 1
            exonAdd = float(int(exonDens) / float(denom))
            intronAdd = float(int(intronDens) / float(denom))
            array = [readLine[0], readLine[1], readLine[2], readLine[3], exonAdd, intronAdd, eiSize, readLine[6]]
            newFile.write(writeFile(array) +'\n')

        readLine = readFile.readline()
        lineNum = lineNum + 1 
        readLine = readLine.strip('\n').split('\t')
                      
        if(readLine[0]!=oldChrom):
            oldChrom = readLine[0]
            if not oldChrom:
                oldChrom = 'done!'
            print "\tNow on chromosome: "+str(oldChrom)
    newFile.close()
    bedFile.close()
    readFile.close()
    return [exonDiscard, intronDiscard]

#REORDER DENSITIES BY CHROMOSOME & GENE
def reorder(output, start):
    shutil.copyfile(str(output)+"\\"+str(start)+"new.txt", str(output)+"\\"+str(start)+"exons.txt")
    oldFile = open(str(output)+"\\"+str(start)+"exons.txt", "r")
    newFile = open(str(output)+"\\"+str(start)+"new.txt", "w")
    geneFile = open(str(output)+"\\"+str(start)+"gene.txt", "r+")
    array = []

    #handles all density information
    line=oldFile.readline()
    while line:
        line = line.split('\t')
        array.append(line)
        line=oldFile.readline()

    array = sorted(array, key=getKey)
    counter = 0
    while counter < len(array):
        newFile.write(writeFile(array[counter]))
        counter = counter + 1
    oldFile.close()
    newFile.close()
    del array

    #handles all genonmic information
    array = []
    line=geneFile.readline()
    while line:
        line = line.split('\t')
        array.append(line)
        line=geneFile.readline()
    geneFile.close()
    geneFile = open(str(output)+"\\"+str(start)+"gene.txt", "w")

    array = sorted(array, key=getKey)
    counter = 0
    while counter < len(array):
        geneFile.write(str(array[counter][0]) +'\t'+ str(array[counter][1]) +'\t'+ str(int(array[counter][2])) +'\t'+ str(int(array[counter][3])) +'\n')
        counter = counter + 1
    geneFile.close()
    del array

#CALCULATE WEIGHTED DENSITIES FOR EACH CHROMOSOME
def total(output, start):
    myFile = open(str(output)+"\\"+str(start)+"new.txt", "r")
    newFile = open(str(output)+"\\"+str(start)+"final.txt", "w")
    discardFile = open(str(output)+"\\"+str(start)+"discard.txt", "a")
    discardedExons = 0
    disGenes = 0

    num_lines = sum(1 for line in open(str(output)+"\\"+str(start)+"new.txt"))
    lineNum = 1
    while lineNum <= num_lines:
        thisLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum).strip('\n').split('\t')
        thisID = thisLine[1]
        basesIntronic = 0
        basesExonic = 0
        weightedE = 0
        weightedI = 0

        newLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum).strip('\n').split('\t')
        while(lineNum<=num_lines and str(newLine[1]) == str(thisID)):
            #newLine[6] == size of E/I
            if(newLine[7] == 'INTRON'):
                basesIntronic = float(basesIntronic) + float(newLine[6])
                weightedI = float(weightedI) + float(float(newLine[5])*float(newLine[6]))
            else:
                basesExonic = float(basesExonic) + float(newLine[6])
                weightedE = float(weightedE) + float(float(newLine[4])*float(newLine[6]))

            lineNum = lineNum + 1
            newLine = linecache.getline(str(output)+"\\"+str(start)+"new.txt", lineNum).strip('\n').split('\t')

        if(basesIntronic != 0 and basesExonic != 0):
            exonDens = (float(weightedE)/float(basesExonic))
            intrDens = (float(weightedI)/float(basesIntronic))
            if(exonDens!=0):
                division = float(float(intrDens)/float(exonDens))
                exonDens = round(exonDens, 3)
                intrDens = round(intrDens, 3)
                division = round(division, 3)
                if(float(division)>=0 and float(division)<=1):
                    newFile.write(str(thisLine[0]) +'\t'+ str(thisLine[1]) +'\t'+ str(thisLine[2]) +'\t'+ str(thisLine[3]) +'\t'+ str(float(exonDens)) +'\t'+ str(float(intrDens)) +'\t'+ str(float(division)) +'\n')
                else:
                    disGenes = disGenes + 1
                    discardFile.write(str(thisLine[0]) +'\t'+ str(thisLine[1]) +'\t'+ str(thisLine[2]) +'\t'+ str(thisLine[3]) +'\t'+ str(float(exonDens)) +'\t'+ str(float(intrDens)) +'\t'+ str(float(division)) +'\tpreMRNA FRAC OUT OF RANGE\n')
            else:
                discardedExons = discardedExons + 1
                discardFile.write(str(thisLine[0]) +'\t'+ str(thisLine[1]) +'\t'+ str(thisLine[2]) +'\t'+ str(thisLine[3]) +'\t'+ str(float(exonDens)) +'\t'+ str(float(intrDens)) +'\tEXONIC DENSITY ZERO\n')

    newFile.close()
    myFile.close()
    return [discardedExons, disGenes]

#ADD SUPPORTING GENOMIC INFORMATION
def getStartStop(output, start):
    gtfFile = open(str(output)+"\\"+str(start)+"gene.txt", "r")
    shutil.copyfile(str(output)+"\\"+str(start)+"final.txt", str(output)+"\\"+str(start)+"exons.txt")
    oldFile = open(str(output)+"\\"+str(start)+"exons.txt", "r")
    newFile = open(str(output)+"\\"+str(start)+"final.txt", "w")

    getLine = oldFile.readline().strip('\n').split('\t')
    gtfLine = gtfFile.readline().strip('\n').split('\t')
    while getLine and gtfLine:
        geneID = getLine[1]
        gtfID = gtfLine[1]

        while str(geneID) != str(gtfID): #until the geneID matches with the line's ID
            gtfLine = gtfFile.readline().strip('\n').split('\t')
            gtfID = gtfLine[1]

	#now write start/stop locations for that gene
        start = gtfLine[2]
        stop = gtfLine[3]
        newFile.write(str(getLine[0]) +'\t'+ str(getLine[1]) +'\t'+ str(getLine[2]) +'\t'+ str(getLine[3]) +'\t'+ str(int(start)) +'\t'+ str(int(stop)) +'\t')
        newFile.write(str(float(getLine[4])) +'\t'+ str(float(getLine[5])) +'\t'+ str(float(getLine[6])) +'\n')
        getLine = oldFile.readline().strip('\n').split('\t')
        if(getLine[0] == ''):
            getLine = False
    gtfFile.close()
    newFile.close()

#CLEANUP TEMPORARY FILES
def cleanup(output, start, finalText): #remove excess files
    totalSum = 0
    outputSum = 0
    
    totalSum = totalSum + os.path.getsize(str(output)+"\\"+str(start)+"exons.txt")
    os.remove(str(output)+"\\"+str(start)+"exons.txt")

    totalSum = totalSum + os.path.getsize(str(output)+"\\"+str(start)+"new.txt")
    os.remove(str(output)+"\\"+str(start)+"new.txt")

    totalSum = totalSum + os.path.getsize(str(output)+"\\"+str(start)+"gene.txt")
    os.remove(str(output)+"\\"+str(start)+"gene.txt")

    totalSum = totalSum + os.path.getsize(str(output)+"\\"+str(start)+"final.txt")
    totalSum = totalSum + os.path.getsize(str(output)+"\\"+str(start)+"discard.txt")
    outputSum = outputSum + os.path.getsize(str(output)+"\\"+str(start)+"final.txt")
    outputSum = outputSum + os.path.getsize(str(output)+"\\"+str(start)+"discard.txt")
    
    if(os.path.isfile(str(output)+"\\"+str(finalText)+".txt")):
        os.remove(str(output)+"\\"+str(finalText)+".txt")
    if(os.path.isfile(str(output)+"\\"+str(finalText)+"_discarded.txt")):
        os.remove(str(output)+"\\"+str(finalText)+"_discarded.txt")

    os.rename(str(output)+"\\"+str(start)+"final.txt", str(output)+"\\"+str(finalText)+".txt")
    os.rename(str(output)+"\\"+str(start)+"discard.txt", str(output)+"\\"+str(finalText)+"_discarded.txt")

    print "Maximum storage: " + str(float(float(totalSum)/float(1048576))) +" MB"
    print "Output storage: " + str(float(float(outputSum)/float(1048576))) +" MB"

#SUPPORTING FUNCTIONS
def listToStr(array):
    return str(array[0]) +'\t'+ str(array[1]) +'\t'+ str(array[2]) +'\t'+ str(array[3]) +'\t'+ str(int(array[4])) +'\t'+ str(int(array[5])) +'\t'+ str(array[6])
def writeFile(array):
    return str(array[0]) +'\t'+ str(array[1]) +'\t'+ str(array[2]) +'\t'+ str(array[3]) +'\t'+ str(float(array[4])) +'\t'+ str(float(array[5])) +'\t'+ str(int(array[6])) +'\t'+ str(array[7])
def getKey(item):
    return [str(item[0]), str(item[1])]

## BEGIN FUNCTION CALLS
print("The current directory of this program is :  "+ str(os.getcwd()))
print""

print("File specifications for BEDGRAPH:")
print("\t-Bedgraph file must be gzipped")
bedGraph = raw_input("Please enter file location/name of the bedgraph file:  ")
if(not os.path.isfile(bedGraph)):
	while(not os.path.isfile(bedGraph)):
		print("Invalid file!")
		bedGraph = raw_input("Please enter file location/name of the bedgraph file:  ")
#get bedGraph name
bedSplit = bedGraph.split('/')
directSize = len(bedSplit)
bedName = bedSplit[-1].split('.')[0]
bedLoc = bedGraph
bedLoc = bedLoc[-1].split('.')
counter = 0
bedLoc = ""
while counter<directSize-1:
    bedLoc+=bedSplit[counter] +'/'
    counter = counter + 1
bedLoc+= str(bedName) + '_chr'
print""

print("File specifications for GTF:")
print("\t-Gtf file must NOT be compressed")
gtfFile = raw_input("Please enter the file location/name of the GTF file: ")
if(not os.path.isfile(gtfFile)):
	while(not os.path.isfile(gtfFile)):
		print("Invalid file!")
		gtfFile = raw_input("Please enter the file location/name of the GTF file: ")
print""

print("Note: This program will be writing several files to reduce memory usage")
output = raw_input("Please enter the output location, use 'sd' for same directory:  ")
if(str(output) == 'sd'):
    output = os.getcwd()
if(not os.path.isdir(output)):
    while(not os.path.isdir(output)):
        print("Invalid directory!")
        output = raw_input("Please enter the output location, use 'sd' for same directory:  ")
        if(str(output) == 'sd'):
            output = os.getcwd()
print""

finalText = raw_input("Please name the final output file:  ")

start = time.time()
print""

print("Retrieving exons...")
getExons(gtfFile, output, start) #extracts exons by genes
print("Exons retrieved!")

print("Removing overlaps...")
loop = True
while loop: #will loop through exons until fully simplified
	loop = condExons(output, start) #condenses overlapping exons by pairs
print("Overlaps removed!")
print (str(int(sum(1 for line in open(str(output)+"\\"+str(start)+"discard.txt")))) +' exons discarded due to overlap')

print("Retrieving introns...")
getIntrons(output, start) #will show intergenic spaces
print("Introns retrieved!")

print("Aligning with bedgraph...")
discarded = densities(bedGraph, output, start) #aligns the genes with bedgraph to find densities for each exon and intron
print("Aligned!")
print str(int(discarded[0])) +" exons not supported by bedgraph file"
print str(int(discarded[1])) +" introns not supported by bedgraph file"

print("Reorganizing data...")
reorder(output, start)
print("Data reorganized!")

print("Calculating densities...")
trashed = total(output, start) #adds up all the individual densities to get the overall density for each gene
print(str(int(trashed[0])) +' genes with zero exonic density')
print(str(int(trashed[1])) +' genes with preMRNA fractions out of range')
print("Densities calculated!")

print("Adding supporting information...")
getStartStop(output, start) #add start/stop locations for each gene
print("Information added!")

print("Cleaning up...")
cleanup(output, start, finalText) #cleans up extra files used
print("Cleaned!")
print""

final = time.time()
print("Please access '"+str(output)+"\\"+str(finalText)+".txt' for output and")
print(str(output)+'\\'+str(finalText)+"_discarded.txt' for discarded information")
print("Time to run:  %s mins" % str(float(float(float(final)-float(start))/float(60))))
print(gtfFile)
raw_input()
