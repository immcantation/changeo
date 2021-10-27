#!/usr/bin/env python3
"""
Converts TSV files into IgPhyML input files
"""

# Info
__author__ = "Kenneth Hoehn"
from changeo import __version__, __date__

# Imports
import os
import random
import subprocess
import multiprocessing as mp
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent
from time import time
from Bio.Seq import Seq
from functools import partial

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.IO import  printLog, printMessage, printWarning, printError, printDebug
from changeo.Defaults import default_format
from changeo.IO import splitName, getDbFields, getFormatOperators, getOutputHandle, getOutputName
from changeo.Alignment import RegionDefinition
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs


def correctMidCodonStart(scodons, qi, debug):
    """
    Find and mask split codons

    Arguments:
      scodons (list): list of codons in IMGT sequence.
      qi (str) : input sequence.
      spos (int) : starting position of IMGT sequence in input sequence.
      debug (bool) : print debugging statements.

    Returns:
      tuple: (modified input sequence, modified starting position of IMGT sequence in input sequence).
    """
    spos = 0   
    for i in range(0, len(scodons)):
        printDebug("%s %s" % (scodons[i], qi[0:3]), debug)
        if scodons[i] != "...":
            if scodons[i][0:2] == "..":
                scodons[i] = "NN" + scodons[i][2]
                #sometimes IMGT will just cut off first letter if non-match, at which point we"ll just want to mask the
                #first codon in the IMGT seq, other times it will be legitimately absent from the query, at which point
                #we have to shift the frame. This attempts to correct for this by looking at the next codon over in the
                #alignment
                if scodons[i][2:3] != qi[2:3] or scodons[i + 1] != qi[3:6]:
                    qi = "NN" + qi
                spos = i
                break
            elif scodons[i][0] == ".":
                scodons[i] = "N" + scodons[i][1:3]
                if scodons[i][1:3] != qi[1:3] or scodons[i+1] != qi[3:6]:
                    qi = "N" + qi
                spos = i
                break
            else:
                spos = i
                break

    return qi, spos


def checkFrameShifts(receptor, oqpos, ospos, log, debug):
    """
    Checks whether a frameshift occured in a sequence

    Arguments:
      receptor (changeo.Receptor.Receptor): Receptor object.
      oqpos (int) : position of interest in input sequence.
      ospos (int) : position of interest in IMGT sequence.
      log (dict) : log of information for each sequence.
      debug (bool) : print debugging statements.
    """   
    frameshifts = 0
    for ins in range(1, 3):
        ros = receptor.sequence_input
        ris = receptor.sequence_imgt
        psite = receptor.v_seq_start - 1 + oqpos*3
        pisite = ospos * 3
        if (psite + 3 + ins) < len(ros) and (pisite + 3) < len(ris):
        #cut out 1 or 2 nucleotides downstream of offending codon
            receptor.sequence_input = ros[0:(psite + 3)] + ros[(psite + 3 + ins):]
            receptor.sequence_imgt = ris[0:(pisite + 3)] + ris[(pisite + 3):]

            # Debug sequence modifications
            printDebug(ros, debug)
            printDebug(receptor.sequence_input, debug)
            printDebug(ris, debug)
            printDebug(receptor.sequence_imgt, debug)
            printDebug("RUNNING %d" % ins, debug)

            mout = maskSplitCodons(receptor, recursive=True)
            if mout[1]["PASS"]:
                #if debug:
                receptor.sequence_input = ros
                receptor.sequence_imgt = ris
                frameshifts += 1
                printDebug("FRAMESHIFT of length %d!" % ins, debug)
                log["FAIL"] = "SINGLE FRAME-SHIFTING INSERTION"
                break
            else:
                receptor.sequence_input = ros
                receptor.sequence_imgt = ris

    return frameshifts


def findAndMask(receptor, scodons, qcodons, spos, s_end, qpos, log, debug, recursive=False):
    """
    Find and mask split codons

    Arguments:
      receptor (changeo.Receptor.Receptor): Receptor object.
      scodons (list): list of codons in IMGT sequence
      qcodons (list): list of codons in input sequence
      spos (int): starting position of IMGT sequence in input sequence
      s_end (int): end of IMGT sequence
      qpos (int): starting position of input sequence in IMGT sequence
      log (dict): log of information for each sequence
      debug (bool): print debugging statements?
      recursive (bool): was this function called recursively?
    """
    frameshifts = 0
    while spos < s_end and qpos < len(qcodons):
        if debug:
            print(scodons[spos] + "\t" + qcodons[qpos])
        if scodons[spos] == "..." and qcodons[qpos] != "...": #if IMGT gap, move forward in imgt
            spos += 1
        elif scodons[spos] == qcodons[qpos]: # if both are the same, move both forward
            spos += 1
            qpos += 1
        elif qcodons[qpos] == "N": # possible that SEQ-IMGT ends on a bunch of Ns
            qpos += 1
            spos += 1
        else: # if not the same, mask IMGT at that site and scan forward until you find a codon that matches next site
            if debug:
                print("checking %s at position %d %d" % (scodons[spos], spos, qpos))
            ospos=spos
            oqpos=qpos
            spos += 1
            qpos += 1
            while spos < s_end and scodons[spos] == "...": #possible next codon is just a gap
                spos += 1
            while qpos < len(qcodons) and spos < s_end and scodons[spos] != qcodons[qpos]:
                printDebug("Checking " + scodons[spos]+ "\t" + qcodons[qpos], debug)
                qpos += 1
            if qcodons[qpos-1] == scodons[ospos]: #if codon in previous position is equal to original codon, it was preserved
                qpos -= 1
                spos = ospos
                printDebug("But codon was apparently preserved", debug)
                if "IN-FRAME" in log:
                    log["IN-FRAME"] = log["IN-FRAME"] + "," +  str(spos)
                else:
                    log["IN-FRAME"] = str(spos)
            elif qpos >= len(qcodons) and spos < s_end:
                printDebug("FAILING MATCH", debug)
                log["PASS"] = False #if no match for the adjacent codon was found, something"s up.
                log["FAIL"] = "FAILED_MATCH_QSTRING:"+str(spos)
                #figure out if this was due to a frame-shift by repeating this method but with an edited input sequence
                if not recursive:
                    frameshifts += checkFrameShifts(receptor, oqpos, ospos, log, debug)

            elif spos >= s_end or qcodons[qpos] != scodons[spos]:
                scodons[ospos] = "NNN"
                if spos >= s_end:
                    printDebug("Masked %s at position %d, at end of subject sequence" % (scodons[ospos], ospos), debug)
                    if "END-MASKED" in log:
                        log["END-MASKED"] = log["END-MASKED"] + "," + str(spos)
                    else:
                        log["END-MASKED"] = str(spos)
                else:
                    printDebug("Masked %s at position %d, but couldn't find upstream match" % (scodons[ospos], ospos), debug)
                    log["PASS"]=False
                    log["FAIL"]="FAILED_MATCH:"+str(spos)
            elif qcodons[qpos] == scodons[spos]:
                printDebug("Masked %s at position %d" % (scodons[ospos], ospos), debug)
                scodons[ospos] = "NNN"
                if "MASKED" in  log:
                    log["MASKED"] = log["MASKED"] + "," + str(spos)
                else:
                    log["MASKED"] = str(spos)
            else:
                log["PASS"] = False
                log["FAIL"] = "UNKNOWN"


def maskSplitCodons(receptor, recursive=False, mask=True):
    """
    Identify junction region by IMGT definition.

    Arguments:
      receptor (changeo.Receptor.Receptor): Receptor object.
      recursive (bool) : was this method part of a recursive call?
      mask (bool) : mask split codons for use with igphyml?

    Returns:
      str: modified IMGT gapped sequence.
      log: dict of sequence information
    """
    debug = False
    qi = receptor.sequence_input
    si = receptor.sequence_imgt
    log = OrderedDict()
    log["ID"]=receptor.sequence_id
    log["CLONE"]=receptor.clone
    log["PASS"] = True
    if debug:
        print(receptor.sequence_id)
    # adjust starting position of query sequence
    qi = qi[(receptor.v_seq_start - 1):]

    #tally where --- gaps are in IMGT sequence and remove them for now
    gaps = []
    ndotgaps = []
    nsi = ""
    for i in range(0,len(si)):
        if si[i] == "-":
            gaps.append(1)
            ndotgaps.append(1)
        else:
            gaps.append(0)
            nsi = nsi + si[i]
            if si[i] != ".":
                ndotgaps.append(0)

    #find any gaps not divisible by three
    curgap = 0
    for i in ndotgaps:
        if i == 1:
            curgap += 1
        elif i == 0 and curgap != 0:
            if curgap % 3 != 0 :
                printDebug("Frame-shifting gap detected! Refusing to include sequence.", debug)
                log["PASS"] = False
                log["FAIL"] = "FRAME-SHIFTING DELETION"
                log["SEQ_IN"] = receptor.sequence_input
                log["SEQ_IMGT"] = receptor.sequence_imgt
                log["SEQ_MASKED"] = receptor.sequence_imgt
                return receptor.sequence_imgt, log
            else:
                curgap = 0
    si = nsi

    scodons = [si[i:i + 3] for i in range(0, len(si), 3)]

    # deal with the fact that it's possible to start mid-codon
    qi,spos = correctMidCodonStart(scodons, qi, debug)

    qcodons = [qi[i:i + 3] for i in range(0, len(qi), 3)]

    frameshifts = 0
    s_end = 0 #adjust for the fact that IMGT sequences can end on gaps
    for i in range(spos, len(scodons)):
        if scodons[i] != "..." and len(scodons[i]) == 3 and scodons[i] != "NNN":
            s_end = i

    printDebug("%i:%i:%s" % (s_end, len(scodons), scodons[s_end]), debug)
    s_end += 1
    qpos = 0

    if mask:
        findAndMask(receptor, scodons, qcodons, spos, s_end, qpos, log, debug, recursive)

    if not log["PASS"] and not recursive:
        log["FRAMESHIFTS"] = frameshifts
    if len(scodons[-1]) != 3:
        if scodons[-1] == ".." or scodons[-1] == ".":
            scodons[-1] = "..."
        else:
            scodons[-1] = "NNN"
        if "END-MASKED" in log:
            log["END-MASKED"] = log["END-MASKED"] + "," + str(len(scodons))
        else:
            log["END-MASKED"] = str(spos)

    concatenated_seq = Seq("")
    for i in scodons:
        concatenated_seq += i

    # add --- gaps back to IMGT sequence
    ncon_seq = ""
    counter = 0
    for i in gaps:
        #print(str(i) + ":" + ncon_seq)
        if i == 1:
            ncon_seq = ncon_seq + "."
        elif i == 0:
            ncon_seq = ncon_seq + concatenated_seq[counter]
            counter += 1
    ncon_seq = ncon_seq + concatenated_seq[counter:]
    concatenated_seq = ncon_seq
    log["SEQ_IN"] = receptor.sequence_input
    log["SEQ_IMGT"] = receptor.sequence_imgt
    log["SEQ_MASKED"] = concatenated_seq

    return concatenated_seq, log


def unAmbigDist(seq1, seq2, fbreak=False):
    """
    Calculate the distance between two sequences counting only A,T,C,Gs

    Arguments:
      seq1 (str): sequence 1
      seq2 (str): sequence 2
      fbreak (bool): break after first difference found?

    Returns:
      int: number of ACGT differences.
    """
    if len(seq1) != len(seq2):
        printError("Sequences are not the same length! %s %s" % (seq1, seq2))

    dist = 0
    for i in range(0,len(seq1)):
        if seq1[i] != "N" and seq1[i] != "-" and seq1[i] != ".":
            if seq2[i] != "N" and seq2[i] != "-" and seq2[i] != ".":
                if seq1[i] != seq2[i]:
                    dist += 1
                    if fbreak:
                        break

    return dist


def deduplicate(useqs, receptors, log=None, meta_data=None, delim=":"):
    """
    Collapses identical sequences

    Argument:
      useqs (dict): unique sequences within a clone. maps sequence to index in Receptor list.
      receptors (dict): receptors within a clone (index is value in useqs dict).
      log (collections.OrderedDict): log of sequence errors.
      meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data.
      meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data.
      delim (str): delimited to use when appending meta_data.

    Returns:
      list: deduplicated receptors within a clone.
    """

    keys = list(useqs.keys())
    join = {} # id -> sequence id to join with (least ambiguous chars)
    joinseqs = {} # id -> useq to join with (least ambiguous chars)
    ambigchar = {} #sequence id -> number ATCG nucleotides
    for i in range(0,len(keys)-1):
        for j in range(i+1,len(keys)):
            ki = keys[i]
            kj = keys[j]
            if meta_data is None:
                ski = keys[i]
                skj = keys[j]
            else:
                ski, cid = keys[i].split(delim)
                skj, cid = keys[j].split(delim)
            ri = receptors[useqs[ki]]
            rj = receptors[useqs[kj]]
            dist = unAmbigDist(ski, skj, True)
            m_match = True
            if meta_data is not None:
                matches = 0
                for m in meta_data:
                    if ri.getField(m) == rj.getField(m) and m != "DUPCOUNT":
                        matches += 1
                m_match = (matches == len(meta_data))
            if dist == 0 and m_match:
                ncounti = ki.count("A") + ki.count("T") + ki.count("G") + ki.count("C")
                ncountj = kj.count("A") + kj.count("T") + kj.count("G") + kj.count("C")
                ambigchar[useqs[ki]] = ncounti
                ambigchar[useqs[kj]] = ncountj
                # this algorithm depends on the fact that all sequences are compared pairwise, and all are zero
                # distance from the sequence they will be collapse to.
                if ncountj > ncounti:
                    nci = 0
                    if useqs[ki] in join:
                        nci = ambigchar[join[useqs[ki]]]
                    if nci < ncountj:
                        join[useqs[ki]] = useqs[kj]
                        joinseqs[ki] = kj
                else:
                    ncj = 0
                    if useqs[kj] in join:
                        ncj = ambigchar[join[useqs[kj]]]
                    if ncj < ncounti:
                        join[useqs[kj]] = useqs[ki]
                        joinseqs[kj] = ki

    # loop through list of joined sequences and collapse
    keys = list(useqs.keys())
    for k in keys:
        if useqs[k] in join:
            rfrom = receptors[useqs[k]]
            rto = receptors[join[useqs[k]]]
            rto.dupcount += rfrom.dupcount
            if log is not None:
                log[rfrom.sequence_id]["PASS"] = False
                log[rfrom.sequence_id]["DUPLICATE"] = True
                log[rfrom.sequence_id]["COLLAPSETO"] = joinseqs[k]
                log[rfrom.sequence_id]["COLLAPSEFROM"] = k
                log[rfrom.sequence_id]["FAIL"] = "Collapsed with " + rto.sequence_id
            del useqs[k]

    return useqs


def hasPTC(sequence):
    """
    Determines whether a PTC exits in a sequence

    Arguments:
      sequence (str): IMGT gapped sequence in frame 1.

    Returns:
      int: negative if not PTCs, position of PTC if found.

    """
    ptcs = ("TAA", "TGA", "TAG", "TRA", "TRG", "TAR", "TGR", "TRR")
    for i in range(0, len(sequence), 3):
            if sequence[i:(i+3)] in ptcs:
                return i
    return -1



def rmCDR3(sequences, clones):
    """
    Remove CDR3 from all sequences and germline of a clone

     Arguments:
       sequences (list): list of sequences in clones.
       clones (list): list of Receptor objects.
    """

    for i in range(0,len(sequences)):
        imgtar = clones[i].getField("imgtpartlabels")
        germline = clones[i].getField("germline_imgt_d_mask")
        nseq = []
        nimgtar = []
        ngermline = []
        ncdr3 = 0
        #print("imgtarlen: " + str(len(imgtar)))
        #print("seqlen: " + str(len(sequences[i])))
        #print("germline: " + str(len(germline)))
        #if len(germline) < len(sequences[i]):
        #    print("\n" + str(clones[i].sequence_id))
        #    print("\n " + str((sequences[i])) )
        #    print("\n" + str((germline)))
        for j in range(0,len(imgtar)):
            if imgtar[j] != 108:
                nseq.append(sequences[i][j])
                if j < len(germline):
                    ngermline.append(germline[j])
                nimgtar.append(imgtar[j])
            else:
                ncdr3 += 1
        clones[i].setField("imgtpartlabels",nimgtar)
        clones[i].setField("germline_imgt_d_mask", "".join(ngermline))
        sequences[i] = "".join(nseq)
        #print("Length: " + str(ncdr3))


def characterizePartitionErrors(sequences, clones, meta_data):
    """
    Characterize potential mismatches between IMGT labels within a clone

    Arguments:
      sequences (list): list of sequences in clones.
      clones (list): list of Receptor objects.
      meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data.
        
    Returns:
      tuple: tuple of length four containing a list of IMGT positions for first sequence in clones,
             the germline sequence of the first receptor in clones, the length of the first sequence in clones,
             and the number of sequences in clones.
    """
    sites = len(sequences[0])
    nseqs = len(sequences)
    imgtar = clones[0].getField("imgtpartlabels")
    germline = clones[0].getField("germline_imgt_d_mask")
    if germline == "":
        germline = clones[0].getField("germline_imgt")

    correctseqs = False

    for seqi in range(0, len(sequences)):
        i = sequences[seqi]
        if len(i) != sites or len(clones[seqi].getField("imgtpartlabels")) != len(imgtar):
            correctseqs = True

    if correctseqs:
        maxlen = sites
        maximgt = len(imgtar)
        for j in range(0,len(sequences)):
            if len(sequences[j]) > maxlen:
                maxlen = len(sequences[j])
            if len(clones[j].getField("imgtpartlabels")) > maximgt:
                imgtar = clones[j].getField("imgtpartlabels")
                maximgt = len(imgtar)
        sites = maxlen
        for j in range(0,len(sequences)):
            cimgt = clones[j].getField("imgtpartlabels")
            seqdiff = maxlen - len(sequences[j])
            imgtdiff = len(imgtar)-len(cimgt)
            sequences[j] = sequences[j] + "N"*(seqdiff)
            last = cimgt[-1]
            cimgt.extend([last]*(imgtdiff))
            clones[j].setField("imgtpartlabels",cimgt)

    if meta_data is not None:
        meta_data_ar = meta_data[0].split(",")
    for c in clones:
        if meta_data is not None:
            c.setField(meta_data[0],c.getField(meta_data_ar[0]))
            for m in range(1,len(meta_data_ar)):
                st = c.getField(meta_data[0])+":"+c.getField(meta_data_ar[m])
                c.setField(meta_data[0],st)

        if len(c.getField("imgtpartlabels")) != len(imgtar):
            printError("IMGT assignments are not the same within clone %d!\n" % c.clone,False)
            printError(c.getField("imgtpartlabels"),False)
            printError("%s\n%d\n" % (imgtar,j),False)
            for j in range(0, len(sequences)):
                printError("%s\n%s\n" % (sequences[j],clones[j].getField("imgtpartlabels")),False)
            printError("ChangeO file needs to be corrected")
        for j in range(0,len(imgtar)):
            if c.getField("imgtpartlabels")[j] != imgtar[j]:
                printError("IMGT assignments are not the same within clone %d!\n" % c.clone, False)
                printError(c.getField("imgtpartlabels"), False)
                printError("%s\n%d\n" % (imgtar, j))

    #Resolve germline if there are differences, e.g. if reconstruction was done before clonal clustering
    resolveglines = False
    for c in clones:
        ngermline = c.getField("germline_imgt_d_mask")
        if ngermline == "":
            ngermline = c.getField("germline_imgt")
        if ngermline != germline:
            resolveglines = True

    if resolveglines:
        printError("%s %s" % ("Predicted germlines are not the same among sequences in the same clone.",
                              "Be sure to cluster sequences into clones first and then predict germlines using --cloned"))

    if sites > (len(germline)):
        seqdiff = sites - len(germline)
        germline = germline + "N" * (seqdiff)

    if sites % 3 != 0:
        printError("number of sites must be divisible by 3! len: %d, clone: %s , id: %s, seq: %s" %(len(sequences[0]),\
        clones[0].clone,clones[0].sequence_id,sequences[0]))
    return imgtar, germline, sites, nseqs

def outputSeqPartFiles(out_dir, useqs_f, meta_data, clones, collapse, nseqs, delim, newgerm, conseqs, duplicate, imgt):
    """
    Create intermediate sequence alignment and partition files for IgPhyML output

    Arguments:
        out_dir (str): directory for sequence files.
        useqs_f (dict): unique sequences mapped to ids.
        meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data.
        clones (list) : list of receptor objects.
        collapse (bool) : deduplicate sequences.
        nseqs (int): number of sequences.
        delim (str) : delimiter for extracting metadata from ID.
        newgerm (str) : modified germline of clonal lineage.
        conseqs (list) : consensus sequences.
        duplicate (bool) : duplicate sequence if only one in a clone.
        imgt (list) : IMGT numbering of clonal positions .
    """

    # bootstrap these data if desired
    lg = len(newgerm)
    sites = range(0, lg)

    transtable = clones[0].sequence_id.maketrans(" ", "_")
    outfile = os.path.join(out_dir, "%s.fasta" % clones[0].clone)
    with open(outfile, "w") as clonef:
        if collapse:
            for seq_f, num in useqs_f.items():
                seq = seq_f
                cid = ""
                if meta_data is not None:
                    seq, cid = seq_f.split(delim)
                    cid = delim + cid.replace(":", "_")
                sid = clones[num].sequence_id.translate(transtable) + cid
                clonef.write(">%s\n%s\n" % (sid.replace(":","-"), seq.replace(".", "-")))
                if len(useqs_f) == 1 and duplicate:
                    if meta_data is not None:
                        if meta_data[0] == "DUPCOUNT":
                            cid = delim + "0"
                    sid = clones[num].sequence_id.translate(transtable) + "_1" + cid
                    clonef.write(">%s\n%s\n" % (sid.replace(":","-"), seq.replace(".", "-")))
        else:
            for j in range(0, nseqs):
                cid = ""
                if meta_data is not None:
                    meta_data_list = []
                    for m in meta_data:
                        meta_data_list.append(clones[j].getField(m).replace(":", "_"))
                    cid = delim + str(delim.join(meta_data_list))
                sid = clones[j].sequence_id.translate(transtable) + cid
                clonef.write(">%s\n%s\n" % (sid.replace(":","-"), conseqs[j].replace(".", "-")))
                if nseqs == 1 and duplicate:
                    if meta_data is not None:
                        if meta_data[0] == "DUPCOUNT":
                            cid = delim + "0"
                    sid = clones[j].sequence_id.translate(transtable)+"_1" + cid
                    clonef.write(">%s\n%s\n" % (sid.replace(":","-"), conseqs[j].replace(".", "-")))

        germ_id = ["GERM"]
        if meta_data is not None:
            for i in range(1,len(meta_data)):
                germ_id.append("GERM")
        clonef.write(">%s_%s\n" % (clones[0].clone,"_".join(germ_id)))
        for i in range(0, len(newgerm)):
            clonef.write("%s" % newgerm[i].replace(".","-"))
        clonef.write("\n")

    #output partition file
    partfile = os.path.join(out_dir, "%s.part.txt" % clones[0].clone)
    with open(partfile, "w") as partf:
        partf.write("%d %d\n" % (2, len(newgerm)))
        partf.write("FWR:IMGT\n")
        partf.write("CDR:IMGT\n")
        partf.write("%s\n" % (clones[0].v_call.split("*")[0]))
        partf.write("%s\n" % (clones[0].j_call.split("*")[0]))
        partf.write(",".join(map(str, imgt)))
        partf.write("\n")


def outputIgPhyML(clones, sequences, meta_data=None, collapse=False, ncdr3=False, logs=None,
                  fail_writer=None, out_dir=None, min_seq=1):
    """
    Create intermediate sequence alignment and partition files for IgPhyML output

    Arguments:
      clones (list): receptor objects within the same clone.
      sequences (list): sequences within the same clone (share indexes with clones parameter).
      meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data
      collapse (bool): if True collapse identical sequences.
      ncdr3 (bool): if True remove CDR3
      logs (dict): contains log information for each sequence
      out_dir (str): directory for output files.
      fail_writer (changeo.IO.TSVWriter): failed sequences writer object.
      min_seq (int): minimum number of data sequences to include.

    Returns:
      int: number of clones.
    """
    s = ""
    delim = "_"
    duplicate = True # duplicate sequences in clones with only 1 sequence?
    imgtar, germline, sites, nseqs = characterizePartitionErrors(sequences, clones, meta_data)

    tallies = []
    for i in range(0, sites, 3):
        tally = 0
        for j in range(0, nseqs):
            if sequences[j][i:(i + 3)] != "...":
                tally += 1
        tallies.append(tally)

    newseqs = []  # remove gap only sites from observed data
    newgerm = []
    imgt = []
    for j in range(0, nseqs):
        for i in range(0, sites, 3):
            if i == 0:
                newseqs.append([])
            if tallies[i//3] > 0:
                newseqs[j].append(sequences[j][i:(i+3)])
    lcodon = ""
    for i in range(0, sites, 3):
        if tallies[i//3] > 0:
            newgerm.append(germline[i:(i+3)])
            lcodon=germline[i:(i+3)]
            imgt.append(imgtar[i])

    if len(lcodon) == 2:
        newgerm[-1] = newgerm[-1] + "N"
    elif len(lcodon) == 1:
        newgerm[-1] = newgerm[-1] + "NN"

    if ncdr3:
        ngerm = []
        nimgt = []
        for i in range(0, len(newseqs)):
            nseq = []
            ncdr3 = 0
            for j in range(0, len(imgt)):
                if imgt[j] != 108:
                    nseq.append(newseqs[i][j])
                    if i == 0:
                        ngerm.append(newgerm[j])
                        nimgt.append(imgt[j])
                else:
                    ncdr3 += 1
            newseqs[i] = nseq
        newgerm = ngerm
        imgt = nimgt
            #print("Length: " + str(ncdr3))

    useqs_f = OrderedDict()
    conseqs = []
    for j in range(0, nseqs):
        conseq = "".join([str(seq_rec) for seq_rec in newseqs[j]])
        if meta_data is not None:
            meta_data_list = []
            for m in range(0,len(meta_data)):
                if isinstance(clones[j].getField(meta_data[m]), str):
                    clones[j].setField(meta_data[m],clones[j].getField(meta_data[m]).replace("_", ""))
                meta_data_list.append(str(clones[j].getField(meta_data[m])))
            conseq_f = "".join([str(seq_rec) for seq_rec in newseqs[j]])+delim+":".join(meta_data_list)
        else:
            conseq_f = conseq
        if conseq_f in useqs_f and collapse:
            clones[useqs_f[conseq_f]].dupcount += clones[j].dupcount
            logs[clones[j].sequence_id]["PASS"] = False
            logs[clones[j].sequence_id]["FAIL"] = "Duplication of " + clones[useqs_f[conseq_f]].sequence_id
            logs[clones[j].sequence_id]["DUPLICATE"]=True
            if fail_writer is not None:
                fail_writer.writeReceptor(clones[j])
        else:
            useqs_f[conseq_f] = j
        conseqs.append(conseq)

    if collapse:
        useqs_f = deduplicate(useqs_f, clones, logs, meta_data, delim)

    if collapse and len(useqs_f) < min_seq:
        for seq_f, num in useqs_f.items():
            logs[clones[num].sequence_id]["FAIL"] = "Clone too small: " + str(len(useqs_f))
            logs[clones[num].sequence_id]["PASS"] = False
        return -len(useqs_f)
    elif not collapse and len(conseqs) < min_seq:
        for j in range(0, nseqs):
            logs[clones[j].sequence_id]["FAIL"] = "Clone too small: " + str(len(conseqs))
            logs[clones[j].sequence_id]["PASS"] = False
        return -len(conseqs)

    # Output fasta file of masked, concatenated sequences
    outputSeqPartFiles(out_dir, useqs_f, meta_data, clones, collapse, nseqs,
                       delim, newgerm, conseqs, duplicate, imgt)

    if collapse:
        return len(useqs_f)
    else:
        return nseqs


def maskCodonsLoop(r, clones, cloneseqs, logs, fails, out_args, fail_writer, mask=True):
    """
    Masks codons split by alignment to IMGT reference

    Arguments:
      r (changeo.Receptor.Receptor): receptor object for a particular sequence.
      clones (list): list of receptors.
      cloneseqs (list): list of masked clone sequences.
      logs (dict): contains log information for each sequence.
      fails (dict): counts of various sequence processing failures.
      out_args (dict): arguments for output preferences.
      fail_writer (changeo.IO.TSVWriter): failed sequences writer object.

    Returns:
      0: returns 0 if an error occurs or masking fails.
      1: returns 1 masking succeeds
    """
    if r.clone is None:
        printError("Cannot export datasets until sequences are clustered into clones.")
    if r.dupcount is None:
        r.dupcount = 1
    fails["rec_count"] += 1
    #printProgress(rec_count, rec_count, 0.05, start_time)
    ptcs = hasPTC(r.sequence_imgt)
    gptcs = hasPTC(r.getField("germline_imgt_d_mask"))
    if gptcs >= 0:
        log = OrderedDict()
        log["ID"] = r.sequence_id
        log["CLONE"] = r.clone
        log["SEQ_IN"] = r.sequence_input
        log["SEQ_IMGT"] = r.sequence_imgt
        logs[r.sequence_id] = log
        logs[r.sequence_id]["PASS"] = False
        logs[r.sequence_id]["FAIL"] = "Germline PTC"
        fails["seq_fail"] += 1
        fails["germlineptc"] += 1
        return 0

    if r.functional and ptcs < 0:
        #If IMGT regions are provided, record their positions
        rd = RegionDefinition(r.junction_length, amino_acid=False)
        regions = rd.getRegions(r.sequence_imgt)
        if regions["cdr3_imgt"] != "" and regions["cdr3_imgt"] is not None:
            simgt = regions["fwr1_imgt"] + regions["cdr1_imgt"] + regions["fwr2_imgt"] + regions["cdr2_imgt"] + \
                    regions["fwr3_imgt"] + regions["cdr3_imgt"] + regions["fwr4_imgt"]
            if len(simgt) < len(r.sequence_imgt):
                r.fwr4_imgt = r.fwr4_imgt +  ("."*(len(r.sequence_imgt) - len(simgt)))
                simgt = regions["fwr1_imgt"] + regions["cdr1_imgt"] + regions["fwr2_imgt"] + \
                        regions["cdr2_imgt"] + regions["fwr3_imgt"] + regions["cdr3_imgt"] + regions["fwr4_imgt"]
            imgtpartlabels = [13]*len(regions["fwr1_imgt"]) + [30]*len(regions["cdr1_imgt"]) + [45]*len(regions["fwr2_imgt"]) + \
                             [60]*len(regions["cdr2_imgt"]) + [80]*len(regions["fwr3_imgt"]) + [108] * len(regions["cdr3_imgt"]) + \
                             [120] * len(regions["fwr4_imgt"])
            r.setField("imgtpartlabels", imgtpartlabels)

            if len(r.getField("imgtpartlabels")) != len(r.sequence_imgt) or simgt != r.sequence_imgt:
                log = OrderedDict()
                log["ID"] = r.sequence_id
                log["CLONE"] = r.clone
                log["SEQ_IN"] = r.sequence_input
                log["SEQ_IMGT"] = r.sequence_imgt
                logs[r.sequence_id] = log
                logs[r.sequence_id]["PASS"] = False
                logs[r.sequence_id]["FAIL"] = "FWR/CDR error"
                logs[r.sequence_id]["FWRCDRSEQ"] = simgt
                fails["seq_fail"] += 1
                fails["region_fail"] += 1
                return 0

        elif regions["fwr3_imgt"] != "" and regions["fwr3_imgt"] is not None:
            simgt = regions["fwr1_imgt"] + regions["cdr1_imgt"] + regions["fwr2_imgt"] + regions["cdr2_imgt"] + \
                    regions["fwr3_imgt"]
            nseq = r.sequence_imgt[len(simgt):len(r.sequence_imgt)]

            if len(simgt) < len(r.sequence_imgt):
                simgt = regions["fwr1_imgt"] + regions["cdr1_imgt"] + regions["fwr2_imgt"] + \
                        regions["cdr2_imgt"] + regions["fwr3_imgt"] + nseq
            imgtpartlabels = [13] * len(regions["fwr1_imgt"]) + [30] * len(regions["cdr1_imgt"]) + [45] * len(
                regions["fwr2_imgt"]) + \
                             [60] * len(regions["cdr2_imgt"]) + [80] * len(regions["fwr3_imgt"]) + \
                             [108] * int(len(nseq))
            r.setField("imgtpartlabels", imgtpartlabels)
            if len(r.getField("imgtpartlabels")) != len(r.sequence_imgt) or simgt != r.sequence_imgt:
                log = OrderedDict()
                log["ID"] = r.sequence_id
                log["CLONE"] = r.clone
                log["SEQ_IN"] = r.sequence_input
                log["SEQ_IMGT"] = r.sequence_imgt
                logs[r.sequence_id] = log
                logs[r.sequence_id]["PASS"] = False
                logs[r.sequence_id]["FAIL"] = "FWR/CDR error"
                logs[r.sequence_id]["FWRCDRSEQ"] = simgt
                fails["seq_fail"] += 1
                fails["region_fail"] += 1
                return 0
        else:
            #imgt_warn = "\n! IMGT FWR/CDR sequence columns not detected.\n! Cannot run CDR/FWR partitioned model on this data.\n"
            imgtpartlabels = [0] * len(r.sequence_imgt)
            r.setField("imgtpartlabels", imgtpartlabels)
        mout = maskSplitCodons(r, mask=mask)
        mask_seq = mout[0]
        ptcs = hasPTC(mask_seq)
        if ptcs >= 0:
            printWarning("Masked sequence suddenly has a PTC.. %s\n" % r.sequence_id)
            mout[1]["PASS"] = False
            mout[1]["FAIL"] = "PTC_ADDED_FROM_MASKING"
        logs[mout[1]["ID"]] = mout[1]
        if mout[1]["PASS"]:
            #passreads += r.dupcount
            if r.clone in clones:
                clones[r.clone].append(r)
                cloneseqs[r.clone].append(mask_seq)
            else:
                clones[r.clone] = [r]
                cloneseqs[r.clone] = [mask_seq]
            return 1
        else:
            if out_args["failed"]:
                fail_writer.writeReceptor(r)
            fails["seq_fail"] += 1
            fails["failreads"] += r.dupcount
            if mout[1]["FAIL"] == "FRAME-SHIFTING DELETION":
                fails["del_fail"] += 1
            elif mout[1]["FAIL"] == "SINGLE FRAME-SHIFTING INSERTION":
                fails["in_fail"] += 1
            else:
                fails["other_fail"] += 1
    else:
        log = OrderedDict()
        log["ID"] = r.sequence_id
        log["CLONE"] = r.clone
        log["PASS"] = False
        log["FAIL"] = "NONFUNCTIONAL/PTC"
        log["SEQ_IN"] = r.sequence_input
        logs[r.sequence_id] = log
        if out_args["failed"]:
            fail_writer.writeReceptor(r)
        fails["seq_fail"] += 1
        fails["nf_fail"] += 1

    return 0

# Run IgPhyML on outputed data
def runIgPhyML(outfile, igphyml_out, clone_dir, nproc=1, optimization="lr", omega="e,e", kappa="e", motifs="FCH",
               hotness="e,e,e,e,e,e", oformat="tab", nohlp=False, asr=-1, clean="none"):
    """
    Run IgPhyML on outputted data

    Arguments:
      outfile (str): Output file name.
      igphymlout (str): igphyml output file
      nproc (int): Number of threads to parallelize IgPhyML across
      optimization (str): Optimize combination of topology (t) branch lengths (l) and parameters (r) in IgPhyML.
      omega (str): omega optimization in IgPhyML (--omega)
      kappa (str): kappa optimization in IgPhyML (-t)
      motifs (str): motifs to use in IgPhyML (--motifs)
      hotness (str): motif in IgPhyML (--hotness)
      oformat (str): output format for IgPhyML (tab or txt)
      nohlp (bool): If True, only estimate GY94 trees and parameters
      clean (str): delete intermediate files? (none, all)
    """
    osplit = outfile.split(".")
    outrep = ".".join(osplit[0:(len(osplit)-1)]) + "_gy.tsv"
    gyout = outfile + "_igphyml_stats_gy.txt"

    gy_args = ["igphyml", "--repfile", outfile, "-m", "GY", "--run_id", "gy", "--outrep", outrep, "--threads",
               str(nproc),"--outname",gyout]


    hlp_args = ["igphyml","--repfile", outrep, "-m", "HLP", "--run_id", "hlp", "--threads", str(nproc), "-o",
                optimization, "--omega", omega, "-t", kappa, "--motifs", motifs, "--hotness", hotness, "--oformat",
                oformat, "--outname", igphyml_out]

    if asr >= 0:
        hlp_args.append("--ASRc")
        hlp_args.append(str(asr))

    log = OrderedDict()
    log["START"] = "IgPhyML GY94 tree estimation"
    printLog(log)

    try: #check for igphyml executable
        subprocess.check_output(["igphyml"])
    except:
        printError("igphyml not found :-/")
    try: #get GY94 starting topologies
        p = subprocess.check_output(gy_args)
    except subprocess.CalledProcessError as e:
        print(" ".join(gy_args))
        print('error>', e.output, '<')
        printError("GY94 tree building in IgPhyML failed")

    log = OrderedDict()
    log["START"] = "IgPhyML HLP analysis"
    log["OPTIMIZE"] = optimization
    log["TS/TV"] = kappa
    log["wFWR,wCDR"] = omega
    log["MOTIFS"] = motifs
    log["HOTNESS"] = hotness
    log["NPROC"] = nproc
    printLog(log)

    if not nohlp:
        try: #estimate HLP parameters/trees
            p = subprocess.check_output(hlp_args)
        except subprocess.CalledProcessError as e:
            print(" ".join(hlp_args))
            print('error>', e.output, '<')
            printError("HLP tree building failed")

    log = OrderedDict()
    log["OUTPUT"] = igphyml_out
    if oformat == "tab":
        igf = open(igphyml_out)
        names = igf.readline().split("\t")
        vals = igf.readline().split("\t")
        for i in range(3,len(names)-1):
            log[names[i]] = round(float(vals[i]),2)
    printLog(log)

    if clean != "none":
        log = OrderedDict()
        log["START"] = "CLEANING"
        log["SCOPE"] = clean
        printLog(log)
        todelete = open(outrep)
        for line in todelete:
            line = line.rstrip("\n")
            line = line.rstrip("\r")
            lsplit = line.split("\t")
            if len(lsplit) == 4:
                os.remove(lsplit[0])
                os.remove(lsplit[1])
                os.remove(lsplit[3])
        todelete.close()
        os.remove(outrep)
        os.remove(outfile)
        os.remove(gyout)
        cilog = outrep + "_igphyml_CIlog.txt_hlp"
        if os.path.isfile(cilog):
            os.remove(cilog)
        if oformat == "tab":
            os.rmdir(clone_dir)
        else:
            printWarning("Using --clean all with --oformat txt will not delete all tree file results.\n"
                         "You'll have to do that yourself.")
        log = OrderedDict()
        log["END"] = "IgPhyML analysis"
        printLog(log)


# Note: Collapse can give misleading dupcount information if some sequences have ambiguous characters at polymorphic sites
def buildTrees(db_file, meta_data=None, target_clones=None, collapse=False, ncdr3=False, nmask=False,
               sample_depth=-1, min_seq=1,append=None, igphyml=False, nproc=1, optimization="lr", omega="e,e",
               kappa="e", motifs="FCH", hotness="e,e,e,e,e,e", oformat="tab", clean="none", nohlp=False,
               asr=-1, format=default_format, out_args=default_out_args):
    """
    Masks codons split by alignment to IMGT reference, then produces input files for IgPhyML

    Arguments:
      db_file (str): input tab-delimited database file.
      meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data
      target_clones (str): List of clone IDs to analyze.
      collapse (bool): if True collapse identical sequences.
      ncdr3 (bool): if True remove all CDR3s.
      nmask (bool): if False, do not attempt to mask split codons
      sample_depth (int): depth of subsampling before deduplication
      min_seq (int): minimum number of sequences per clone
      append (str): column name to append to sequence_id
      igphyml (bool): If True, run IgPhyML on outputted data
      nproc (int) : Number of threads to parallelize IgPhyML across
      optimization (str): Optimize combination of topology (t) branch lengths (l) and parameters (r) in IgPhyML.
      omega (str): omega optimization in IgPhyML (--omega)
      kappa (str): kappa optimization in IgPhyML (-t)
      motifs (str): motifs to use in IgPhyML (--motifs)
      hotness (str): motif in IgPhyML (--hotness)
      oformat (str): output format for IgPhyML (tab or txt)
      clean (str): delete intermediate files? (none, all)
      nohlp (bool): If True, only estimate GY94 trees and parameters
      format (str): input and output format.
      out_args (dict): arguments for output preferences.

    Returns:
      dict: dictionary of output pass and fail files.
    """
    # Print parameter info
    log = OrderedDict()
    log["START"] = "BuildTrees"
    log["FILE"] = os.path.basename(db_file)
    log["COLLAPSE"] = collapse
    printLog(log)

    # Open output files
    out_label = "lineages"

    pass_handle = getOutputHandle(db_file,
                                  out_label=out_label,
                                  out_dir=out_args["out_dir"],
                                  out_name= out_args["out_name"],
                                  out_type="tsv")

    igphyml_out = None
    if igphyml:
        igphyml_out = getOutputName(db_file, out_label="igphyml-pass",
                                    out_dir=out_args["out_dir"],
                                    out_name=out_args["out_name"],
                                    out_type=oformat)

    dir_name, __ = os.path.split(pass_handle.name)

    if out_args["out_name"] is None:
        __, clone_name, __ = splitName(db_file)
    else:
        clone_name = out_args["out_name"]
    if dir_name is None:
        clone_dir = clone_name
    else:
        clone_dir = os.path.join(dir_name, clone_name)
    if not os.path.exists(clone_dir):
        os.makedirs(clone_dir)

    # Format options
    try:
        reader, writer, __ = getFormatOperators(format)
    except ValueError:
        printError("Invalid format %s." % format)
    out_fields = getDbFields(db_file, reader=reader)

    # open input file
    handle = open(db_file, "r")
    records = reader(handle)

    fail_handle, fail_writer = None, None
    if out_args["failed"]:
        fail_handle = getOutputHandle(db_file,
                                      out_label="lineages-fail",
                                      out_dir=out_args["out_dir"],
                                      out_name=out_args["out_name"],
                                      out_type=out_args["out_type"])
        fail_writer = writer(fail_handle, fields=out_fields)

    cloneseqs = {}
    clones = {}
    logs = OrderedDict()
    fails = {"rec_count":0, "seq_fail":0, "nf_fail":0, "del_fail":0, "in_fail":0, "minseq_fail":0,
             "other_fail":0, "region_fail":0, "germlineptc":0, "fdcount":0, "totalreads":0,
             "passreads":0, "failreads":0}

    # Mask codons split by indels
    start_time = time()
    printMessage("Correcting frames and indels of sequences", start_time=start_time, width=50)

    #subsampling loop
    init_clone_sizes = {}
    big_enough = []
    all_records = []
    found_no_funct = False
    for r in records:
        if r.functional is None:
            r.functional = True
            if found_no_funct is False:
                printWarning("FUNCTIONAL column not found.")
                found_no_funct = True
        all_records.append(r)
        if r.clone in init_clone_sizes:
            init_clone_sizes[r.clone] += 1
        else:
            init_clone_sizes[r.clone] = 1

    for r in all_records:
        if target_clones is None or r.clone in target_clones:
            if init_clone_sizes[r.clone] >= min_seq:
               big_enough.append(r)

    fails["totalreads"] = len(all_records)
    #fails["minseq_fail"] =  len(all_records) - len(big_enough)         

    if len(big_enough) == 0:
        printError("\n\nNo sequences found that match specified criteria.",1)

    if sample_depth > 0:
        random.shuffle(big_enough)

    total = 0
    for r in big_enough:
        if r.functional is None:
            r.functional = True
            if found_no_funct is False:
                printWarning("FUNCTIONAL column not found.")
                found_no_funct = True

        r.sequence_id = r.sequence_id.replace(",","-") #remove commas from sequence ID
        r.sequence_id = r.sequence_id.replace(":","-") #remove colons from sequence ID
        r.sequence_id = r.sequence_id.replace(",","-") #remove commas from sequence ID
        r.sequence_id = r.sequence_id.replace(")","-") #remove parenthesis from sequence ID
        r.sequence_id = r.sequence_id.replace("(","-") #remove parenthesis from sequence ID
        if(meta_data is not None):
            for m in range(0,len(meta_data)):
                md = r.getField(meta_data[m])
                md = md.replace(",","-") #remove commas from metadata
                md = md.replace(":","-") #remove colons from metadata
                md = md.replace(",","-") #remove commas from metadata
                md = md.replace(")","-") #remove parenthesis from metadata
                md = md.replace("(","-") #remove parenthesis from metadata
                r.setField(meta_data[m],md)

        if append is not None:
            if append is not None:
                for m in append:
                    r.sequence_id = r.sequence_id + "_" + r.getField(m)
        total += maskCodonsLoop(r, clones, cloneseqs, logs, fails, out_args, fail_writer,
                                mask = not nmask)
        if total == sample_depth:
            break

    # Start processing clones
    clonesizes = {}
    pass_count, nclones = 0, 0
    printMessage("Processing clones", start_time=start_time, width=50)
    for k in clones.keys():
        if len(clones[str(k)]) < min_seq:
            for j in range(0, len(clones[str(k)])):
                logs[clones[str(k)][j].sequence_id]["FAIL"] = "Clone too small: " + str(len(cloneseqs[str(k)]))
                logs[clones[str(k)][j].sequence_id]["PASS"] = False
            clonesizes[str(k)] = -len(cloneseqs[str(k)])
        else:
            clonesizes[str(k)] = outputIgPhyML(clones[str(k)], cloneseqs[str(k)], meta_data=meta_data, collapse=collapse,
                                           ncdr3=ncdr3, logs=logs, fail_writer=fail_writer,
                                           out_dir=clone_dir, min_seq=min_seq)

        #If clone is too small, size is returned as a negative
        if clonesizes[str(k)] > 0:
            nclones += 1
            pass_count += clonesizes[str(k)]
        else:
            fails["seq_fail"] -= clonesizes[str(k)]
            fails["minseq_fail"]  -= clonesizes[str(k)]
    fail_count = fails["rec_count"] - pass_count

    # End clone processing
    printMessage("Done", start_time=start_time, end=True, width=50)

    log_handle = None
    if out_args["log_file"] is not None:
        log_handle = open(out_args["log_file"], "w")
        for j in logs.keys():
            printLog(logs[j], handle=log_handle)

    pass_handle.write(str(nclones)+"\n")
    for key in sorted(clonesizes, key=clonesizes.get, reverse=True):
        #print(key + "\t" + str(clonesizes[key]))
        outfile = os.path.join(clone_dir, "%s.fasta" % key)
        partfile = os.path.join(clone_dir, "%s.part.txt" % key)
        if clonesizes[key] > 0:
            germ_id = ["GERM"]
            if meta_data is not None:
                for i in range(1, len(meta_data)):
                    germ_id.append("GERM")
            pass_handle.write("%s\t%s\t%s_%s\t%s\n" % (outfile, "N", key,"_".join(germ_id), partfile))

    handle.close()
    output = {"pass": None, "fail": None}
    if pass_handle is not None:
        output["pass"] = pass_handle.name
        pass_handle.close()
    if fail_handle is not None:
        output["fail"] = fail_handle.name
        fail_handle.close()
    if log_handle is not None:
        log_handle.close()

    #printProgress(rec_count, rec_count, 0.05, start_time)
    log = OrderedDict()
    log["OUTPUT"] = os.path.basename(pass_handle.name) if pass_handle is not None else None
    log["RECORDS"] = fails["totalreads"]
    log["INITIAL_FILTER"] = fails["rec_count"]
    log["PASS"] = pass_count
    log["FAIL"] = fail_count
    log["NONFUNCTIONAL"] = fails["nf_fail"]
    log["FRAMESHIFT_DEL"] = fails["del_fail"]
    log["FRAMESHIFT_INS"] = fails["in_fail"]
    log["CLONETOOSMALL"] = fails["minseq_fail"]
    log["CDRFWR_ERROR"] = fails["region_fail"]
    log["GERMLINE_PTC"] = fails["germlineptc"]
    log["OTHER_FAIL"] = fails["other_fail"]

    if collapse:
        log["DUPLICATE"] = fail_count - fails["seq_fail"]
    log["END"] = "BuildTrees"
    printLog(log)

    #Run IgPhyML on outputted data?
    if igphyml:
        runIgPhyML(pass_handle.name, igphyml_out=igphyml_out, clone_dir=clone_dir, nproc=nproc,
                   optimization=optimization, omega=omega, kappa=kappa, motifs=motifs,
                   hotness=hotness, oformat=oformat, nohlp=nohlp,clean=clean,asr=asr)

    return output


def getArgParser():
    """
    Defines the ArgumentParser

    Returns:
        argparse.ArgumentParser: argument parsers.
    """
    # Define input and output field help message
    fields = dedent(
             """
             output files:
                 <folder>
                   folder containing fasta and partition files for each clone.
                 lineages
                    successfully processed records.
                 lineages-fail
                    database records failed processing.
                 igphyml-pass
                    parameter estimates and lineage trees from running IgPhyML, if specified

             required fields:
                 sequence_id, sequence, sequence_alignment,
                 germline_alignment_d_mask or germline_alignment,
                 v_call, j_call, clone_id, v_sequence_start
              """)

    # Parent parser
    parser_parent = getCommonArgParser(out_file=False, log=True, format=True)

    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[parser_parent],
                            formatter_class=CommonHelpFormatter, add_help=False)

    group = parser.add_argument_group("sequence processing arguments")
    group.add_argument("--collapse", action="store_true", dest="collapse",
                        help="""If specified, collapse identical sequences before exporting to fasta.""")
    group.add_argument("--ncdr3", action="store_true", dest="ncdr3",
                        help="""If specified, remove CDR3 from all sequences.""")
    group.add_argument("--nmask", action="store_true", dest="nmask",
                       help="""If specified, do not attempt to mask split codons.""")
    group.add_argument("--md", nargs="+", action="store", dest="meta_data",
                       help="""List of fields to containing metadata to include in output fasta file 
                            sequence headers.""")
    group.add_argument("--clones", nargs="+", action="store", dest="target_clones",
                       help="""List of clone IDs to output, if specified.""")
    group.add_argument("--minseq", action="store", dest="min_seq", type=int, default=1,
                       help="""Minimum number of data sequences. Any clones with fewer than the specified
                            number of sequences will be excluded.""")
    group.add_argument("--sample", action="store", dest="sample_depth", type=int, default=-1,
                       help="""Depth of reads to be subsampled (before deduplication).""")
    group.add_argument("--append", nargs="+", action="store", dest="append",
                       help="""List of columns to append to sequence ID to ensure uniqueness.""")

    igphyml_group = parser.add_argument_group("IgPhyML arguments (see igphyml -h for details)")
    igphyml_group.add_argument("--igphyml", action="store_true", dest="igphyml",
                               help="""Run IgPhyML on output?""")
    igphyml_group.add_argument("--nproc", action="store", dest="nproc", type=int, default=1,
                               help="""Number of threads to parallelize IgPhyML across.""")
    igphyml_group.add_argument("--clean", action="store", choices=("none", "all"),
                               dest="clean", type=str, default="none",
                               help="""Delete intermediate files? 
                               none: leave all intermediate files; all: delete all intermediate files.""")
    igphyml_group.add_argument("--optimize", action="store", dest="optimization", type=str, default="lr",
                               choices=("n","r","l","lr","tl","tlr"),
                               help="""Optimize combination of topology (t) branch lengths (l) and parameters (r), or 
                               nothing (n), for IgPhyML.""")
    igphyml_group.add_argument("--omega", action="store", dest="omega", type=str, default="e,e",
                               help="""Omega parameters to estimate for FWR,CDR respectively: 
                               e = estimate, ce = estimate + confidence interval, or numeric value""")
    igphyml_group.add_argument("-t", action="store", dest="kappa", type=str, default="e",
                               help="""Kappa parameters to estimate: 
                               e = estimate, ce = estimate + confidence interval, or numeric value""")
    igphyml_group.add_argument("--motifs", action="store", dest="motifs", type=str,
                               default="WRC_2:0,GYW_0:1,WA_1:2,TW_0:3,SYC_2:4,GRS_0:5",
                               help="""Which motifs to estimate mutability.""")
    igphyml_group.add_argument("--hotness", action="store", dest="hotness", type=str, default="e,e,e,e,e,e",
                               help="""Mutability parameters to estimate: 
                               e = estimate, ce = estimate + confidence interval, or numeric value""")
    igphyml_group.add_argument("--oformat", action="store", dest="oformat", type=str, default="tab",
                               choices=("tab", "txt"),
                               help="""IgPhyML output format.""")
    igphyml_group.add_argument("--nohlp", action="store_true", dest="nohlp",
                               help="""Don't run HLP model?""")
    igphyml_group.add_argument("--asr", action="store", dest="asr", type=float, default=-1,
                               help="""Ancestral sequence reconstruction interval (0-1).""")

    return parser


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    # Parse command line arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    del args_dict["db_files"]

    # Call main for each input file
    for f in args.__dict__["db_files"]:
        args_dict["db_file"] = f
        buildTrees(**args_dict)