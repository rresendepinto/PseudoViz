from csv import list_dialects
from math import trunc
import os, sys, time, re, subprocess, shutil
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter, length_hint
import argparse
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool   
from contextlib import closing
import time
import pandas as pd
import json


def upper_str(nucl):
    return nucl.upper()

def replace(word, pattern1, subst):
    return(word.replace(pattern1, subst))

def averagepairwiseaaidentity(results_aa):
        identities = []
        avgidentity = 0
        for i in range(1, len(results_aa)):
            identity = 0
            for j in range(len(results_aa[i].seq)):
                if results_aa[i].seq[j] == results_aa[0].seq[j]:
                    identity += 1
            identity = identity / len(results_aa[i])
            identities.append(identity)
        for i in identities:
            avgidentity += i

        avgidentity = avgidentity / len(identities) * 100
        avgidentity = round(avgidentity, 2)
        return str(avgidentity)

def aaidenticalsites(results_aa):
        t = len(results_aa[0].seq)
        
        for i in range(len(results_aa[0].seq)):
            for j in range(1, len(results_aa)):
                if results_aa[0].seq[i] != results_aa[j].seq[i]:
                    t -= 1
                    break
        return str(t)

def alignmentlength(results_aa):
        return str(len(results_aa[0].seq) * 3)

def nosequences(results_aa):
    return str(len(results_aa))

def intersect_strings(cds, utr_exon):
    #works only for first exon
    clean_exon = utr_exon

    while clean_exon != cds[0:len(clean_exon)]:

        clean_exon = clean_exon[1:]

    if clean_exon != cds[0:len(clean_exon)]:

        return None

    else:

        return clean_exon

def intersect_strings_intronless(cds, utr_exon):

    clean_exon = utr_exon

    while cds != clean_exon[0:len(cds)]:

        clean_exon = clean_exon[1:]

    clean_exon = clean_exon[0:len(cds)]

    if clean_exon != cds[0:len(clean_exon)]:

        return None

    else:

        return clean_exon

class Mutation():

    def __init__(self, start, end, type_mut):   

        self.start = start
        self.end = end
        self.type_mut = type_mut

    def get_mut_length(self):

        return(self.end + 1 - self.start)

def find_stop_codon(seq):
    stopcodons = {'TAA', 'TAG', 'TGA'}
    stop_presence=False
    stop_list = []
    for p in range(0, len(seq), 3):
        if p != 0 and p != len(seq) - 3:
            if seq[p : p + 3] in stopcodons:
                stop_presence = True
                stop_list.append(p + 1)

    return(stop_list)

def find_stop_codon2(seq1, refseq2):

    stopcodons = {'TAA', 'TAG', 'TGA'}
    
    stop_list = []
    #find where reference sequence actually starts, i.e. after any possible leading gaps
    start_of_ref = 0
    while refseq2[start_of_ref] == '-':
        start_of_ref += 1  
    #get indices of gaps in the reference sequence
    indices = [i for i, x in enumerate(refseq2) if x == "-" and i > start_of_ref]
    
    #now search the target sequence with gaps for stop codons
    for p in range(0, len(seq1), 3):

        if p != len(seq1) - 3:

            if seq1[p : p + 3] in stopcodons:   
                #add the stop codons to the list but subtract the number of gaps in the reference sequence before the stop codon
                stop_list.append(p + 1 - len([x for x in indices if x < p + 1])) 

    return(stop_list)
                       

def get_pos_exons(list_exs_cds, cds):
    num_exons = len([e for e in list_exs_cds if e.id != 'CDS'])
    pos_exons = {}
    last_end = 0

    for e in list_exs_cds:  # Aligns each exon against the genomic region
        if e.id != 'CDS':

            if e.description.endswith('Exon_1') and e.description.endswith('Exon_' + str(num_exons)):
                #for intronless genes
                exon_seq = intersect_strings_intronless(cds, e.seq)

                if len(exon_seq) == 0:
                    return None
                else:
                    pos_exons[e.id] = (1, len(exon_seq))
                    last_end = len(exon_seq)

            elif e.description.endswith('Exon_1'):

                exon_seq = intersect_strings(cds, e.seq)

                if len(exon_seq) == 0:
                    return None
                else:
                    pos_exons[e.id] = (1, len(exon_seq))
                    last_end = len(exon_seq)

            elif e.description.endswith('Exon_' + str(num_exons)):
                #this excludes 3' utrs 
                pos_exons[e.id] = (last_end + 1, min(len(cds), last_end + len(e.seq)))

            else:

                pos_exons[e.id] = (last_end + 1, last_end + len(e.seq))
                last_end = last_end + len(e.seq)

    return(pos_exons)                       
        
def convert_exon_pos2global_pos(list_exs_cds, cds, exon, exon_pos):

    exon_pos_dict = get_pos_exons(list_exs_cds, cds)

    exon = 'Exon_' + str(exon)

    global_pos = exon_pos + exon_pos_dict[exon][0] -1 

    return(global_pos)

def get_exon_from_globalpos(list_exs_cds, global_pos):

    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)

    exon_pos_dict = get_pos_exons(list_exs_cds, cds)

    exon_determined = 0

    for e in exon_pos_dict.keys():
        
        if global_pos >= exon_pos_dict[e][0] and global_pos <= exon_pos_dict[e][1]:
            
            exon_determined = e 

    if exon_determined != 0:

        return exon_determined






def check_exon_mutations(ref, target):

    mutations = []
    
    #gap opened
    gap = None

    #where is the gap
    gap_which = ''

    for i in range(len(ref)):
        #if difference between sequences
        
        if ref[i] != target[i]:
            
            if gap == None:

                if ref[i] == '-':

                    gap=[i, i]

                    gap_which = 'Reference'

                elif target[i] == '-':

                    gap=[i, i]

                    gap_which = 'Target'

            else:

                if gap_which == 'Reference':

                    if ref[i] == '-':
                        gap[1] = i


                    else:

                        x = Mutation(gap[0], gap[1], 'insertion')
                        mutations.append(x)                        

                        gap=None

                        gap_which = ''

                        if target[i] == '-':
                            
                            gap=[i, i]

                            gap_which = 'Target'

                elif gap_which == 'Target':

                    
                    if target[i] == '-':
                        gap[1] = i


                    else:

                        x = Mutation(gap[0], gap[1], 'deletion')
                        mutations.append(x)

                        gap=None

                        gap_which = ''

                        if ref[i] == '-':
                            
                            gap=[i, i]

                            gap_which = 'Reference'


        elif gap:

            if gap_which == 'Reference':

                x = Mutation(gap[0] + 1, gap[1] + 1, 'insertion')
                mutations.append(x)
                gap = None

            elif gap_which == 'Target':

                x = Mutation(gap[0] + 1, gap[1] + 1, 'deletion')
                mutations.append(x)
                gap = None

        

    
    return mutations

def create_mutations_table(exon_alns, id):
    #this function takes a genomic target class and returns a table with annotated mutation 
    #it iterates over the exon_alns attribute of this class

    ###TODO need to change this to take an exon alns or something like that###
    exon_dict = {}
    
    for exon in exon_alns.keys():
        
        record_ref = SeqRecord(Seq(exon_alns[exon][1].upper()),  id="Reference_species")
        
        record_target = SeqRecord(Seq(exon_alns[exon][0].upper()),  id=id)

        exon_mutations = check_exon_mutations(record_ref.seq, record_target.seq)
        
        exon_dict[exon] = exon_mutations

    list_frameshifts = []   

    for exon in exon_dict.keys():

        if exon_dict[exon]:

            for mut in exon_dict[exon]:

                
        
                list_frameshifts.append([id,exon,mut.type_mut, mut.get_mut_length(), mut.start, mut.end])

    list_frameshifts = [x for x in list_frameshifts if x[3] % 3 != 0]
    table_frameshifts = pd.DataFrame(list_frameshifts, columns=['id', 'exon', 'type','len', 'start', 'end'])
    return(table_frameshifts)

def reverse_modulo(n, x):

    mod = n%x
    rev_mod = x-mod

    return(rev_mod)
        
def create_pstop_table(exon_alns , list_exs_cds, id, table_frameshifts):
    #this functions creates a table with stop codons in the job's main path
    
    t_aligned_exons = [x for x in exon_alns.keys()]
    
    #to find stop codons first we need a sequence with the exons joined but with gaps representing the missing exons
    #so the lack of an exon doesnt alter the reading frame like a frameshift would
    
    exonscat_withgaps = ''
    reference_exonscat_withgaps = ''

    exon_pos_dict = get_pos_exons(list_exs_cds=list_exs_cds, cds=str(list_exs_cds[len(list_exs_cds) - 1].seq))
    
    for e in list_exs_cds:
        
        if e.id != 'CDS':
            
            exon_n = e.description.split('_')[-1]
           
            if exon_n in t_aligned_exons :
                exonscat_withgaps += exon_alns[exon_n][0]
                reference_exonscat_withgaps += exon_alns[exon_n][1]
            else:
                exonscat_withgaps += '-' * (exon_pos_dict[e.id][1] - exon_pos_dict[e.id][0] + 1)
                reference_exonscat_withgaps += '*' * (exon_pos_dict[e.id][1] - exon_pos_dict[e.id][0] + 1)

    
    #essential for correct search of stop codons after frameshift insertions
    if table_frameshifts.empty == False:
        for index, row in table_frameshifts.iterrows():
            if (row[3] % 3 != 0):
                if row[2] == 'insertion':
                    cds=str(list_exs_cds[len(list_exs_cds) - 1].seq)
                    global_pos_frameshift = convert_exon_pos2global_pos(list_exs_cds, cds, row[1], row[4])
                    #like MACSE, it will add a gap to the insertion to keep the sequence in frame
                    exonscat_withgaps = exonscat_withgaps[0:global_pos_frameshift] + reverse_modulo(row[3], 3) * '-' + exonscat_withgaps[global_pos_frameshift:]
                    #to keep in frame the reference sequence also needs to be modified
                    reference_exonscat_withgaps = reference_exonscat_withgaps[0:global_pos_frameshift] + reverse_modulo(row[3], 3) * '-' + reference_exonscat_withgaps[global_pos_frameshift:]
        
    

    #returns a list of positions containing premature stop codons
    pstops = find_stop_codon2(exonscat_withgaps, reference_exonscat_withgaps)
    

    
    list_pstops = []

    

    for pstop in pstops:

        #pstop = pstop - reference_exonscat_withgaps[0:pstop].count('-')

        pstop_exon = get_exon_from_globalpos(list_exs_cds, pstop) 
    
        list_pstops.append([id, pstop_exon, pstop])


    if len(list_pstops) > 0:
        table_pstops = pd.DataFrame(list_pstops, columns=['id', 'exon', 'pos'])
        
    
        return (table_pstops, (exonscat_withgaps, reference_exonscat_withgaps))

    else:

        table_pstops = pd.DataFrame()
        
        return (table_pstops, exonscat_withgaps)

def get_amino_from_cds_dict(cds_dict):
    #not currently used because I cant figure out a way to present aminoacid sequences in exons with frameshifts
    amino_dict = {}

    for sp in cds_dict.keys():

        amino_dict[sp] = (Seq(cds_dict[sp][0]).translate, Seq(cds_dict[sp][0]).translate())

    return amino_dict

def shifted_codons(table_mut, list_exs_cds):
    #outputs the percentage of the cds that is shifted    

    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)
    
    #how shifted the frame is
    frameshift_score = 0

    #length of shifted frame
    total_shifted = 0

    #save position of previous mutation
    previous_mut = 0



    for index, row in table_mut.iterrows():
       
        #get position of mutation in cds instead of in exon
        mut_pos = convert_exon_pos2global_pos(list_exs_cds, cds, row['exon'], row['start'])
        
        #if frame is shifted add the length since previous mutation or since start of cds
        if frameshift_score % 3 != 0:

            total_shifted += mut_pos - previous_mut

        #change the score according to shiftedness of frame
        if row['type'] == "insertion":

            frameshift_score += row['len']

        elif row['type'] == "deletion":

            frameshift_score -= row['len']

        previous_mut = mut_pos
        
        
    #at the end add the lenght of shifted cds since last mutation (if the frame is still shifted)
    if frameshift_score % 3 != 0:
        
        total_shifted += len(cds) - previous_mut
    
    return total_shifted/len(cds) * 100

def truncated_codons(table_pstop, list_exs_cds):
    
    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)

    first_pstop = table_pstop['pos'].min()

    return (len(cds) - first_pstop ) / len(cds) * 100

def absent_exons(exon_alns, list_exs_cds):
    #determines the percentage of cds lost by missing exons
    list_exs = [e for e in list_exs_cds if e.id != 'CDS']
    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)
    
    present_exons = [x for x in exon_alns.keys()]
    
    exon_pos_dict = get_pos_exons(list_exs, cds)

    total_present = 0

    if len(present_exons) != 0 and len(list_exs) == 1:
        return 0

    for e in present_exons:

        #if len(list_exs) > 1:
        total_present += exon_pos_dict['Exon_' + str(e)][1] - exon_pos_dict['Exon_' + str(e)][0] +1
        



    return (1 - total_present/len(cds)) * 100

def alt_pseudoindex(shifted_codons, truncated_codons, absent_exons):
    #pseudoindex created from the first step of the software , the alignment of the exons
    #since the original pseudoindex is created from the MACSE alignment
    pseudoindex = 0

    #first lets check the percentage of shifted codons from frameshifts
    if shifted_codons < 10:
        pseudoindex = max(0, pseudoindex)
    elif shifted_codons > 10 and shifted_codons < 15:
        pseudoindex = max(1, pseudoindex)
    elif shifted_codons > 15 and shifted_codons < 20:
        pseudoindex = max(2, pseudoindex)
    elif shifted_codons > 20 and shifted_codons < 25:
        pseudoindex = max(3, pseudoindex)
    elif shifted_codons > 25 and shifted_codons < 30:
        pseudoindex = max(4, pseudoindex)
    elif shifted_codons > 30:
        pseudoindex = 5
    
    #now lets check the percentage of absent cds from missing exons
    if absent_exons < 10:
        pseudoindex = max(0, pseudoindex)
    elif absent_exons > 10 and absent_exons < 15:
        pseudoindex = max(1, pseudoindex)
    elif absent_exons > 15 and absent_exons < 20:
        absent_exons = max(2, pseudoindex)
    elif absent_exons > 20 and absent_exons < 25:
        pseudoindex = max(3, pseudoindex)
    elif absent_exons > 25 and absent_exons < 30:
        pseudoindex = max(4, pseudoindex)
    elif absent_exons > 30:
        pseudoindex = 5


    #now lets check the percentage of truncated cds
    if truncated_codons < 10:
        pseudoindex = max(0, pseudoindex)
    elif truncated_codons > 10 and truncated_codons < 15:
        pseudoindex = max(1, pseudoindex)
    elif truncated_codons > 15 and truncated_codons < 20:
        pseudoindex = max(2, pseudoindex)
    elif truncated_codons > 20 and truncated_codons < 25:
        pseudoindex = max(3, pseudoindex)
    elif truncated_codons > 25 and truncated_codons < 30:
        pseudoindex = max(4, pseudoindex)
    elif truncated_codons > 30:
        pseudoindex = 5

    return pseudoindex

def pseudoindex_stats(table_frameshifts, table_pstops, exon_alns, list_exs_cds):

        shifted_codons_percent = shifted_codons(table_frameshifts, list_exs_cds)

        if table_pstops.empty:
            truncated_codons_percent = 0
            
         
        else:
            truncated_codons_percent = truncated_codons(table_pstops, list_exs_cds)
            
        
        missing_percent = absent_exons(exon_alns, list_exs_cds)
        

        pseudoindex =  alt_pseudoindex(shifted_codons_percent, truncated_codons_percent, missing_percent)
        
        list_stats = [pseudoindex, shifted_codons_percent, truncated_codons_percent, missing_percent]

        #table_pseudoindex = pd.DataFrame(list_stats, columns=['pseudoindex', 'shifted_frame%', 'truncated_frame%', 'missing_frame%'])

        return list_stats


def general_stats(exon_data, list_exs_cds):

    '''Coding sequence prediction statistics'''

    def avgidentity(exon_alns):
        '''Returns the average alignment identity'''
        a = 0
        for i in exon_alns:
            a += exon_alns[i][6]
        avgidentity = a / len(exon_alns)
        avgidentity = str(round(avgidentity, 2))
        return avgidentity

    def avgexonsize(exon_alns):
        '''Returns the average predicted exon size'''
        b = 0
        for i in exon_alns:
            b += exon_alns[i][8]-exon_alns[i][7] + 1
        avgexonsize = b / len(exon_alns)
        avgexonsize = str(round(avgexonsize, 2))
        return avgexonsize

    def minexon(exon_alns):
        '''Returns a string with the smallest exon and the corresponding size'''
        exsize = []
        for i in exon_alns:
            exsize.append([exon_alns[i][8]-exon_alns[i][7] + 1, i])
        exsize = sorted(exsize)
        return str('Exon ' + str(exsize[0][1]) + ' (' + str(exsize[0][0]) + ' bp)')

    def maxexon(exon_alns):
        '''Returns a string with the largest exon and the corresponding size'''
        exsize = []
        for i in exon_alns:
            exsize.append([exon_alns[i][8]-exon_alns[i][7] + 1, i])
        exsize = sorted(exsize)
        return str('Exon ' + str(exsize[len(exsize) - 1][1]) + ' (' + str(exsize[len(exsize) - 1][0]) + ' bp)')

    def splicesiteintegrity(exon_alns):
        '''Returns a string with the splice site integrity (No. of functional splicing sites/total splicing sites)'''
        sp = 0
        sptotal = 0
        if numtotalexons != 1:  # Intronless genes are not considered
            for i in exon_alns:
                
                if int(i) == 1:
                    if exon_alns[i][4].upper() == 'GT' or exon_alns[i][4].upper() == 'GC':
                        sp += 1
                    sptotal += 1
                elif int(i) == numtotalexons:
                    if exon_alns[i][3].upper() == 'AG':
                        sp += 1
                    sptotal += 1
                else:
                    if exon_alns[i][3].upper() == 'AG':
                        sp += 1
                    if exon_alns[i][4].upper() == 'GT' or exon_alns[i][4].upper() == 'GC':
                        sp += 1
                    sptotal += 2
        
        return str(sp) + '/' + str(sptotal)

    general_stats_dict = {}
    numtotalexons = len(list_exs_cds) - 1
    for target in exon_data.keys():


        try:
            general_stats_dict[target] = [avgidentity(exon_data[target]), avgexonsize(exon_data[target]), minexon(exon_data[target]), maxexon(exon_data[target]), splicesiteintegrity(exon_data[target]), ','.join([x for x in exon_data[target].keys()])]
        except:
            pass

    return general_stats_dict

def create_pstops_list(df, list_exs_cds, cds):

    list_pstops = []

    exon_pos_dict = get_pos_exons(list_exs_cds, cds)

    for i, row in df.iterrows():
        
        exon_n = get_exon_from_globalpos(list_exs_cds, row['pos'])
        
        exon_pos = row["pos"] - exon_pos_dict[exon_n][0] + 1
        list_pstops.append(f'{row["id"]}_{exon_n}_{exon_pos}')
        list_pstops.append(f'{row["id"]}_{exon_n}_{exon_pos + 1}')
        list_pstops.append(f'{row["id"]}_{exon_n}_{exon_pos + 2}')

    return list_pstops

def translate_seq(seq):
    seq_nt = Seq(seq)
    seq_aa = seq_nt.translate()
    return str(seq_aa)

