import os
import csv
import re

#Prints directory of this script for convenience
print (os.getcwd()) 
#Sets working directory
os.chdir(r"/Users/file_path/file_directory") 

#Opens input csv file located in working directory
with open('input.csv','r') as csv_file: 
    lines = csv_file.readlines()

#Appends imported csv values to the lists below
PeptideSequences = []
PSI = []
for line in lines:
    line = line.strip('\n')
    data = line.split(',')
    PeptideSequences.append(data[0]) #Peptide sequence from first column of csv
    PSI.append(data[1]) #Protein Stability Index from second column of csv

#Removes column headers imported from csv
PeptideSequences.remove('Peptide Sequence')
PSI.remove('PSI')
print(PeptideSequences[-1]) #check that list of sequence properly imports by printing last sequence in list

#Defines dictionary of lists used to record amino acid fraction composition of each peptide sequence
composition = {"A":[], "C":[], "D":[], "E":[], "F":[], "G":[], "H":[], "I":[], "K":[], "L":[], 
               "M":[], "N":[], "P":[], "Q":[], "R":[], "S":[], "T":[], "V":[], "W":[], "Y":[]}
#Lists for recording sequence features
number_unique_residues = []
composition_fraction_square_sum = []

#Takes string as input, appends composition of string to "composition" dictionary of lists
#Calculates amino acid fractional composition, number of unique residues, and fraction squared sum simultaneously
def SeqToComp(sequence): 
    unique_amino_acids = 0
    fraction_square = 0
    for key in composition:
        fraction = (sequence.count(key))/len(sequence)
        composition[key].append(fraction)
        if fraction > 0:
            unique_amino_acids += 1
            fraction_square += fraction**2
    number_unique_residues.append(unique_amino_acids)
    composition_fraction_square_sum.append(fraction_square)
    return composition

#Takes string and appends 1 if it contains minimal BAG6 motif
BAG6 = []
def BAG6finder(i): 
    if re.search("[L,I,V,M][L,I,V,M][L,I,V,M][L,I,V,M]", i) or re.search("[L,I,V,M].[L,I,V,M][L,I,V,M][L,I,V,M]", i) or \
    re.search("[L,I,V,M][L,I,V,M].[L,I,V,M][L,I,V,M]", i) or re.search("[L,I,V,M][L,I,V,M][L,I,V,M].[L,I,V,M]", i) or \
    re.search("[L,I,V,M]..[L,I,V,M][L,I,V,M][L,I,V,M]", i) or re.search("[L,I,V,M].[L,I,V,M].[L,I,V,M][L,I,V,M]", i) or \
    re.search("[L,I,V,M].[L,I,V,M][L,I,V,M].[L,I,V,M]", i) or re.search("[L,I,V,M][L,I,V,M]..[L,I,V,M][L,I,V,M]", i) or \
    re.search("[L,I,V,M][L,I,V,M].[L,I,V,M].[L,I,V,M]", i) or re.search("[L,I,V,M][L,I,V,M][L,I,V,M]..[L,I,V,M]", i):
        BAG6.append('1')
    else:
        BAG6.append('0')

#Takes string and appends 1 if it contains minimal amphipathic motif
Amphipathic = []
def Amphipathicfinder(i):
    if re.search("[L,I,V,F,Y,W]..[L,I,V,F,Y,W]...[L,I,V,F,Y,W]", i) or \
    re.search("[L,I,V,F,Y,W]...[L,I,V,F,Y,W]..[L,I,V,F,Y,W]", i):
        Amphipathic.append('1')
    else:
        Amphipathic.append('0')

#Finds repeats of disorder promoting residues and calculates 
#a score equal to the sum of k*(# of k-mer) from k = 2 to k = sequence length
Repeat_score = []
def Repeatfinder(i):
    score = 0
    unique_repeat_count = 0
    repeat_length = len(i)
    for n in range(repeat_length):
        if repeat_length >= 2:
            search = '([A,R,G,Q,E,K,P,S])\\1{' + str(repeat_length-1) + ',}' #'Disorder promoting" amino acids enriched in IDRs
            kmer_count = len(re.findall(search, i))
            score += repeat_length*(kmer_count - unique_repeat_count)
            repeat_length -= 1
            unique_repeat_count += kmer_count - unique_repeat_count
    Repeat_score.append(score)

#Takes string and calculates largest motif length and largest motif score
Motif_max_length = []
Motif_max_score = []
#Pearson correlation coefficients of each amino acid to protein stability index, used in a weighing function
correlation = {"A":0.0758, "C":-0.1190, "D":0.2234, "E":0.3529, "F":-0.4060, "G":0.1434, "H":-0.0237, "I":-0.3295, "K":0.0866, "L":-0.4229, 
               "M":-0.1451, "N":0.0084, "P":0.3057, "Q":0.1639, "R":0.0184, "S":0.2260, "T":0.0687, "V":-0.2118, "W":-0.2719, "Y":-0.2861}
def Motif_max(i):
    current_position = 0
    current_best_length = 0
    current_best_score = 0
    scoreboard = []
    total_positive = 0
    for residue in i:
        #only residues with negative correlation with PSI
        if residue in ('C','H','F','I','L','M','V','W','Y'):
            total_positive += 1
        scoreboard.append(-1*correlation[residue])
    for a in range(total_positive): #iterates window length
        for b in range(len(i)-a): #iterates window position
            if sum(scoreboard[current_position+b:current_position+b+a+1]) > current_best_score:
                current_best_score = sum(scoreboard[current_position+b:current_position+b+a+1])
                current_best_length = a + 1

    Motif_max_length.append(current_best_length)
    Motif_max_score.append(current_best_score)

#Iterates over each sequence in "PeptideSequences" and uses the above functions to calculate features
for i in PeptideSequences:
    SeqToComp(i)
    BAG6finder(i)
    Amphipathicfinder(i)
    Repeatfinder(i)
    Motif_max(i)

#Add extra features to composition for export
composition['BAG6 Motif'] = BAG6
composition['Amphipathic'] = Amphipathic
composition['Fraction Squared Sum'] = composition_fraction_square_sum
composition['Number of Unique Residues'] = number_unique_residues
composition['Repeat Score'] = Repeat_score
composition['Largest Motif Length'] = Motif_max_length
composition['Largest Motif Score'] = Motif_max_score
composition['PSI'] = PSI

#Outputs "composition" dictionary of lists to csv in same folder as script
with open("GNMT Fragments Output.csv", "w", newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(composition.keys())
    writer.writerows(zip(*composition.values()))
