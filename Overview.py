#!/usr/bin/python
import sys
import math
import re

###Program to give quick overview of input DNA, RNA, or protein sequence. 

##Takes both uppercase and lowercase inputs
##Stop codons in AA sequence depicted as " * " 

#DNA or RNA inputs will produce the other along with the protein sequnce. DNA composition given to 2 decimal places.
#A truncated protein sequnce may be produced starting from the first Start codon to the first Stop codon.
#A Molecular weight of both proteins will be provided to 2 decimals of the kDa

#---------------------------------------------------------------------------------
#AA Letter coding: 
#https://meme-suite.org/meme/doc/alphabets.html

#AA Weights key derived from: Amino Acid Abbreviations and Molecular Weights
#https://www.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/#:~:text=The%20average%20molecular%20weight%20of,of%2064%2C000%20grams%20per%20mole.

#AA Codon key derived from: Standard DNA codon table
#https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
#---------------------------------------------------------------------------------

if (len(sys.argv) < 3):
    print("\n Requires 2 arguments: \'-p\' for protein, \'-r\' for RNA, and \'-d\' for DNA. Followed by the nucleotide or amino acid string. \n")
    sys.exit()

count=0
for input in sys.argv:
    if (input == '-d'):
        inputtype = "d"
        typeloc = count
    if (input == '-r'):
        inputtype = "r"
        typeloc = count
    if (input == '-p'):
        inputtype = "p"
        typeloc = count
    count+=1
try: inputtype
except NameError: inputtype = None   
if (inputtype == None):
    sys.exit("\n Error: No acceptable input type given. \n")


##Variables
AA = {
  "TTT": "F",
  "TTC": "F",
  "TTA": "L",
  "TTG": "L",

  "CTT": "L",
  "CTC": "L",
  "CTA": "L",
  "CTG": "L",

  "ATT": "I",
  "ATC": "I",
  "ATA": "I",
  "ATG": "M",

  "GTT": "V",
  "GTC": "V",
  "GTA": "V",
  "GTG": "V",
  
  "TCT": "S",
  "TCC": "S",
  "TCA": "S",
  "TCG": "S",

  "CCT": "P",
  "CCC": "P",
  "CCA": "P",
  "CCG": "P",

  "ACT": "T",
  "ACC": "T",
  "ACA": "T",
  "ACG": "T",

  "GCT": "A",
  "GCC": "A",
  "GCA": "A",
  "GCG": "A",

  "TAT": "Y",
  "TAC": "Y",
  "TAA": "*",
  "TAG": "*",

  "CAT": "H",
  "CAC": "H",
  "CAA": "Q",
  "CAG": "Q",

  "AAT": "N",
  "AAC": "N",
  "AAA": "K",
  "AAG": "K",

  "GAT": "D",
  "GAC": "D",
  "GAA": "E",
  "GAG": "E",

  "TGT": "C",
  "TGC": "C",
  "TGA": "*",
  "TGG": "W",

  "CGT": "R",
  "CGC": "R",
  "CGA": "R",
  "CGG": "R",

  "AGT": "S",
  "AGC": "S",
  "AGA": "R",
  "AGG": "R",

  "GGT": "G",
  "GGC": "G",
  "GGA": "G",
  "GGG": "G",
}
AAweight = {
  "A": 89,
  "R": 174,
  "N": 132,
  "D": 133,
  "B": 133,
  "C": 121,
  "Q": 146,
  "E": 147,
  "Z": 147,
  "G": 75,
  "H": 155,
  "I": 131,
  "L": 131,
  "K": 146,
  "M": 149,
  "F": 165,
  "P": 115,
  "S": 105,
  "T": 119,
  "W": 204,
  "Y": 181,
  "V": 117,
  "X": 110,
  "*": 0
}
MW = 0
Stops = 0
Acont = 0
Tcont = 0
Gcont = 0
Ccont = 0
TruncActive = 0
TruncProtein = ""


##DNA input
if (inputtype == 'd'):
    print("\n DNA input: " + str(sys.argv[typeloc+1]).upper())
    RNAseq = ""
    for nucleotide in sys.argv[typeloc+1]:
        if (nucleotide.upper() == 'A'):
            RNAseq = (RNAseq + "U")
            Acont+=1
            continue
        if (nucleotide.upper() == 'T'):
            RNAseq = (RNAseq + "A")
            Tcont+=1
            continue
        if (nucleotide.upper() == 'G'):
            RNAseq = (RNAseq + "C")
            Gcont+=1
            continue
        if (nucleotide.upper() == 'C'):
            RNAseq = (RNAseq + "G")
            Ccont+=1
            continue
        sys.exit("Error DNA input invalid \n")
    DNAseq = str(sys.argv[typeloc+1])
    Acont=Acont/len(DNAseq)
    Tcont=Tcont/len(DNAseq)
    Gcont=Gcont/len(DNAseq)
    Ccont=Ccont/len(DNAseq)
    print(" A: %.2f  T: %.2f  G: %.2f  C: %.2f \n" % (Acont, Tcont, Gcont, Ccont) )
    print((" RNA Sequence: " + RNAseq + "\n"))
    codons = []
    for group in range(math.floor(len(sys.argv[typeloc+1])/3)):
        codon = sys.argv[typeloc+1][(group*3)] + sys.argv[typeloc+1][(group*3)+1] + sys.argv[typeloc+1][(group*3)+2] 
        codons.append(str(AA[codon.upper()]))
        if (AA[codon.upper()] == "*"):
            Stops+=1
    Proteinseq = ''.join(codons)
    print(" Protein Sequence: " + str(Proteinseq))
        

##RNA input
if (inputtype == 'r'):
    print("\n RNA input: " + str(sys.argv[typeloc+1]).upper() + "\n")
    DNAseq = ""
    for nucleotide in sys.argv[typeloc+1]:
        if (nucleotide.upper() == 'U'):
            DNAseq = (DNAseq + "A")
            Acont+=1
            continue
        if (nucleotide.upper() == 'A'):
            DNAseq = (DNAseq + "T")
            Tcont+=1
            continue
        if (nucleotide.upper() == 'C'):
            DNAseq = (DNAseq + "G")
            Gcont+=1
            continue
        if (nucleotide.upper() == 'G'):
            DNAseq = (DNAseq + "C")
            Ccont+=1
            continue
        sys.exit("Error RNA input invalid \n")
    RNAseq = str(sys.argv[typeloc+1]).upper()
    print((" DNA Sequence: " + DNAseq))
    Acont=Acont/len(DNAseq)
    Tcont=Tcont/len(DNAseq)
    Gcont=Gcont/len(DNAseq)
    Ccont=Ccont/len(DNAseq)
    print(" A: %.2f  T: %.2f  G: %.2f  C: %.2f \n" % (Acont, Tcont, Gcont, Ccont) )
    codons = []
    for group in range(math.floor(len(sys.argv[typeloc+1])/3)):
        codon = DNAseq[(group*3)] + DNAseq[(group*3)+1] + DNAseq[(group*3)+2] 
        codons.append(str(AA[codon.upper()]))
        if (AA[codon.upper()] == "*"):
            Stops+=1
    Proteinseq = ''.join(codons)
    print(" Protein Sequence: " + str(Proteinseq))


##Protein input
if (inputtype == 'p'):
    Proteinseq = sys.argv[typeloc+1].upper()
    print("\n Protein input: " + str(Proteinseq))


##All Sequences continue here with protein info
for AA in Proteinseq:
    MW = MW + AAweight[AA]
MW = MW - (18.015 * (len(Proteinseq)-1))
MW = MW / 1000
print(" Protein length: " + str(len(Proteinseq)-Stops))
print(" Appoximate weight: %.02f KDa \n" % MW)

#Generate Truncated protein
for AA in Proteinseq:
    if (AA == 'M'):
        TruncActive = 1
    if (TruncActive == 1):
        TruncProtein = TruncProtein + AA
        if (AA == '*'):
            TruncActive = 0
            break
if (len(TruncProtein) > 0):
    print(" Truncated Protein: " + TruncProtein)
    print(" Tuncated Protein length: " + str(len(TruncProtein)-1))
    MW = 0
    for AA in TruncProtein:
        MW = MW + AAweight[AA]
    MW = MW - (18.015 * (len(TruncProtein)-2))
    MW = MW / 1000
    print(" Appoximate weight: %.02f KDa \n" % MW)

