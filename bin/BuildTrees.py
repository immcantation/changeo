#!/usr/bin/env python3
"""
Converts TSV files into IgPhyML input files
"""

# Imports
import os
import sys
from argparse import ArgumentParser
from textwrap import dedent
from collections import OrderedDict
from time import time

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.IO import  printLog, getOutputHandle, printProgress
from changeo.Defaults import default_format
from changeo.IO import AIRRReader, ChangeoReader, AIRRWriter, ChangeoWriter, \
                       splitFileName, getDbFields, getFormatOperators
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs

from Bio.Seq import Seq

def maskSplitCodons(receptor):
    """
    Identify junction region by IMGT definition.

    Arguments:
      receptor : Receptor object.

    Returns:
      str : modified IMGT gapped sequence.
    """

    qi = receptor.sequence_input
    si = receptor.sequence_imgt
    log = OrderedDict()
    log['ID']=receptor.sequence_id
    log['CLONE']=receptor.clone
    log['PASS'] = True
   # print(receptor.sequence_id)
    # adjust starting position of query sequence
    qi = qi[(receptor.v_seq_start - 1):]

    # deal with the fact that it's possible to start mid-codon
    scodons = [si[i:i + 3] for i in range(0, len(si), 3)]
    for i in range(0, len(scodons)):
        #print(scodons[i],qi[0:3])
        if scodons[i] != '...':
            if scodons[i][0:2] == '..':
                scodons[i] = "NN"+scodons[i][2]
                #sometimes IMGT will just cut off first letter if non-match, at which point we'll just want to mask the
                #first codon in the IMGT seq, other times it will be legitimately absent from the query, at which point
                #we have to shift the frame. This attempts to correct for this by looking at the next codon over in the
                #alignment
                if scodons[i][2:3] != qi[2:3] or scodons[i + 1] != qi[3:6]:
                    qi = "NN" + qi
                spos = i
                break
            elif scodons[i][0] == '.':
                scodons[i] = "N" + scodons[i][1:3]
                if scodons[i][1:3] != qi[1:3] or scodons[i+1] != qi[3:6]:
                    #print("here!\t"+scodons[i][1:3] +"\t"+ qi[0:2]+":"+scodons[i+1]+"\t"+qi[3:6])
                    qi = "N" + qi
                spos = i
                break
            else:
                spos = i
                break

    qcodons = [qi[i:i + 3] for i in range(0, len(qi), 3)]


    qpos = 0
    for i in range(spos, len(scodons)):
        if scodons[i] != '...':
            qpos += 1

    qpos = 0
    # TODO: for loop with zip()
    while spos < len(scodons) and qpos < len(qcodons):
        #print(scodons[spos] + "\t" + qcodons[qpos])
        if scodons[spos] == '...' and qcodons[qpos] != '...': #if IMGT gap, move forward in imgt
            spos += 1
        elif scodons[spos] == qcodons[qpos]: # if both are the same, move both forward
            spos += 1
            qpos += 1
        else: # if not the same, mask IMGT at that site and scan foward until you find a codon that matches next site
            #print("checking %s at position %d" % (scodons[spos], spos))
            ospos=spos
            spos += 1
            qpos += 1
            while qpos < len(qcodons) and spos < len(scodons) and scodons[spos] != qcodons[qpos]:
                qpos += 1
            if qcodons[qpos-1] == scodons[ospos]: #if codon in previous position is equal to original codon, it was preserved
                qpos -= 1
                spos = ospos
               # print("But codon was apparently preserved")
                #log[str(spos)]="IN-FRAME"
                if 'IN-FRAME' in  log:
                    log['IN-FRAME'] = log['IN-FRAME'] + "," +  str(spos)
                else:
                    log['IN-FRAME'] = str(spos)
            elif qpos >= len(qcodons) and spos < len(scodons):
                log['PASS'] = False
                log['FAIL'] = "FAILED_MATCH_QSTRING:"+str(spos)
            elif spos >= len(scodons) or qcodons[qpos] != scodons[spos]:
                scodons[ospos] = "NNN"
                if spos >= len(scodons):
                #    print("Masked %s at position %d, at end of subject sequence" % (scodons[ospos], ospos))
                    #log[str(spos)] = "END"
                    if 'END-MASKED' in log:
                        log['END-MASKED'] = log['END-MASKED'] + "," + str(spos)
                    else:
                        log['END-MASKED'] = str(spos)
                else:
                 #   print("Masked %s at position %d, but couldn't find upstream match" % (scodons[ospos], ospos))
                    #log[str(spos)] = "FAILED_MATCH"
                    log['PASS']=False
                    log['FAIL']="FAILED_MATCH:"+str(spos)
                    #exit(1)
            elif qcodons[qpos] == scodons[spos]:
               # print("Masked %s at position %d" % (scodons[ospos], ospos))
                scodons[ospos] = "NNN"
                if 'MASKED' in  log:
                    log['MASKED'] = log['MASKED'] + "," + str(spos)
                else:
                    log['MASKED'] = str(spos)

            else:
                log['PASS'] = False
                log['FAIL'] = "UNKNOWN"

    #if scodons[-1] == "." or scodons[-1] == ".." or scodons[-1] == "...":
    #    scodons[-1] = "NNN"
    #    log[str(len(scodons))] = "MASKED"
    if len(scodons[-1]) != 3:
        scodons[-1] = "NNN"
        #log[str(len(scodons))] = "MASKED"
        if 'MASKED' in log:
            log['MASKED'] = log['MASKED'] + "," + str(len(scodons))
        else:
            log['MASKED'] = str(spos)

    concatenated_seq = Seq("")
    for i in scodons:
        concatenated_seq += i

    log['SEQ_IN'] = receptor.sequence_input
    log['SEQ_IMGT'] = receptor.sequence_imgt
    log["SEQ_MASKED"] = concatenated_seq

    return concatenated_seq, log

def outputIgPhyML(clones, sequences, meta_data=None, collapse=False, logs=None, fail_writer=None, out_dir=None):
    """
    Create intermediate sequence alignment and partition files for IgPhyML output

    Arguments:
        clones (list): receptor objects within the same clone.
        sequences (list): sequences within the same clone (share indexes with clones parameter).
        meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data
        collapse (bool): if True collapse identical sequences.
        logs (dict of OrderedDicts): contains log information for each sequence
        out_dir (str): directory for output files.
        fail_writer (handle): file of failed sequences
    Returns:
        int: number of clones.
    """
    duplicate = True # duplicate sequences in clones with only 1 sequence?
    sites = len(sequences[0])
    s=""
    for i in sequences:
        if len(i) != sites:
            print("Sequences within clone %s are not the same length!" % clones[0].clone)
            for s in sequences:
                print(s)
            exit(1)

    nseqs = len(sequences)
    germline = clones[0].getField("germline_imgt_d_mask")

    if len(sequences[0]) % 3 != 0:
        print("number of sites must be divisible by 3! len: %d, clone: %s , seq: %s" %(len(sequences),clones[0].clone,sequences[0]))
        exit(1)
    tallies = []
    for i in range(0, sites, 3):
        tally = 0
        for j in range(0, nseqs):
            if sequences[j][i:(i+3)] != "...":
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
    lcodon = ''
    for i in range(0, sites, 3):
        if tallies[i//3] > 0:
            newgerm.append(germline[i:(i+3)])
            lcodon=germline[i:(i+3)]
            imgt.append(i//3)

    if len(lcodon) == 2:
        newgerm[-1]=newgerm[-1]+"N"
    elif len(lcodon) == 1:
        newgerm[-1] = newgerm[-1] + "NN"

    transtable = clones[0].sequence_id.maketrans(" ","_")

    delim = "_"
    #useqs = {}
    useqs_f = {}
    conseqs = []
    for j in range(0, nseqs):
        conseq = "".join([str(seq_rec) for seq_rec in newseqs[j]])
        if meta_data is not None:
            #print("Meta data:",meta_data[0],meta_data)
            #print(clones[j].getField(meta_data[0]))
            conseq_f = "".join([str(seq_rec) for seq_rec in newseqs[j]])+delim+str(clones[j].getField(meta_data[0]))
        else:
            conseq_f = conseq
        if conseq_f in useqs_f and collapse:
            logs[clones[j].sequence_id]['PASS'] = False
            logs[clones[j].sequence_id]['FAIL'] = "Duplication of " + clones[useqs_f[conseq_f]].sequence_id
            logs[clones[j].sequence_id]['DUPLICATE']=True
            if fail_writer is not None:
                fail_writer.writeReceptor(clones[j])
        else:
            #useqs[conseq] = j
            useqs_f[conseq_f] = j
        conseqs.append(conseq)

    # Output fasta file of masked, concatenated sequences
    outfile = out_dir + "/" + clones[0].clone + ".fa"
    clonef = open(outfile, 'w')
    if collapse:
        for seq_f, num in useqs_f.items():
            seq = seq_f
            cid = ""
            if meta_data is not None:
                seq, cid = seq_f.split(delim)
                cid = delim + cid
            sid = clones[num].sequence_id.translate(transtable) + cid
            print(">%s\n%s" % (sid, seq), file=clonef)
            if len(useqs_f) == 1 and duplicate:
                sid = clones[num].sequence_id.translate(transtable) + "_1" + cid
                print(">%s\n%s" % (sid, seq), file=clonef)
    else:
        for j in range(0, nseqs):
            cid = ""
            if meta_data is not None:
                cid = delim+str(clones[j].getField(meta_data[0]))
            sid = clones[j].sequence_id.translate(transtable)+cid
            print(">%s\n%s" % (sid, conseqs[j]), file=clonef)
            if nseqs == 1 and duplicate:
                sid = clones[j].sequence_id.translate(transtable)+"_1"+cid
                print(">%s\n%s" % (sid, conseqs[j]), file=clonef)

    print(">%s_GERM" % clones[0].clone, file=clonef)
    for i in range(0, len(newgerm)):
        print("%s" % newgerm[i], end='', file=clonef)
    print("\n", end='', file=clonef)
    clonef.close()

    #output partition file
    partfile = out_dir + "/" + clones[0].clone + ".part.txt"
    partf = open(partfile, 'w')
    print("%d %d" % (2, len(newgerm)), file=partf)
    print("FWR:IMGT", file=partf)
    print("CDR:IMGT", file=partf)
    print("%s" % (clones[0].v_call.split("*")[0]),file=partf)
    print("%s" % (clones[0].j_call.split("*")[0]),file=partf)
    print(",".join(map(str, imgt)), file=partf)

    if collapse:
        return len(useqs_f)
    else:
        return nseqs


def buildTrees(db_file, meta_data=None, collapse=False, format=default_format, out_args=default_out_args):
    """
    Masks codons split by alignment to IMGT reference, then produces input files for IgPhyML

    Arguments:
        db_file (str): input tab-delimited database file.
        meta_data (str): Field to append to sequence IDs. Splits identical sequences with different meta_data
        collapse (bool): if True collapse identical sequences.
        format (str): input and output format.
        out_args (dict): arguments for output preferences.
    Returns:
        None: outputs files for IgPhyML.

    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'BuildTrees'
    log['FILE'] = os.path.basename(db_file)
    log['COLLAPSE'] = collapse
    printLog(log)

    start_time = time()
    pass_handle = getOutputHandle(db_file,
                                  out_label='lineages',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type='tsv')

    dir_name, __ = os.path.split(pass_handle.name)

    if out_args['out_name'] is None:
        clone_name, __ = splitFileName(db_file)
    else:
        clone_name = out_args['out_name']
    # clone_dir = outdir/out_name
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
        sys.exit('Error:  Invalid format %s' % format)
    out_fields = getDbFields(db_file, reader=reader)

    # open input file
    handle = open(db_file, 'r')
    records = reader(handle)

    fail_handle, fail_writer = None, None
    if out_args['failed']:
        fail_handle = getOutputHandle(db_file,
                                      out_label='lineages-fail',
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_args['out_type'])
        fail_writer = writer(fail_handle, fields=out_fields)

    cloneseqs = {}
    clones = {}
    logs = OrderedDict()
    rec_count = 0
    seq_fail = 0
    for r in records:
        rec_count += 1
        #printProgress(rec_count, rec_count, 0.05, start_time)
        if r.functional:
            mout = maskSplitCodons(r)
            mask_seq = mout[0]
            logs[mout[1]['ID']]=mout[1]
            if mout[1]['PASS']:
                if r.clone in clones:
                    clones[r.clone].append(r)
                    cloneseqs[r.clone].append(mask_seq)
                else:
                    clones[r.clone] = [r]
                    cloneseqs[r.clone] = [mask_seq]
            else:
                if out_args['failed']:
                    fail_writer.writeReceptor(r)
                seq_fail += 1

        else:
            log = OrderedDict()
            log['ID'] = r.sequence_id
            log['CLONE'] = r.clone
            log['PASS'] = False
            log['FAIL'] = 'NONFUNCTIONAL'
            log['SEQ_IN'] = r.sequence_input
            logs[r.sequence_id]=log
            if out_args['failed']:
                fail_writer.writeReceptor(r)
            seq_fail += 1

    clonesizes = {}
    pass_count = 0
    for k in clones.keys():
        clonesizes[str(k)] = outputIgPhyML(clones[str(k)], cloneseqs[str(k)], meta_data=meta_data, collapse=collapse,
                                           logs=logs,fail_writer=fail_writer,out_dir=clone_dir)
        pass_count += clonesizes[str(k)]
    fail_count = rec_count - pass_count

    log_handle = None
    if out_args['log_file'] is not None:
        log_handle = open(out_args['log_file'], 'w')
        for j in logs.keys():
            printLog(logs[j], handle=log_handle)

    # TODO: changeo console log
    print(len(clonesizes), file=pass_handle)
    for key in sorted(clonesizes, key=clonesizes.get, reverse=True):
        #print(key + "\t" + str(clonesizes[key]))
        outfile = clone_dir + "/" + key + ".fa"
        partfile = clone_dir + "/" + key + ".part.txt"
        pass_handle.write("%s\t%s\t%s\t%s\n" % (outfile, "N", key+"_GERM", partfile))

    handle.close()
    output = {'pass': None, 'fail': None}
    if pass_handle is not None:
        output['pass'] = pass_handle.name
        pass_handle.close()
    if fail_handle is not None:
        output['fail'] = fail_handle.name
        fail_handle.close()
    if log_handle is not None:
        log_handle.close()

    #printProgress(rec_count, rec_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name) if pass_handle is not None else None
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['MASKFAIL'] = seq_fail
    if collapse:
        log['DUPLICATE'] = fail_count - seq_fail
    log['END'] = 'BuildTrees'
    printLog(log)

def getArgParser():
    """
    Defines the ArgumentParser

    Returns:
        argparse.ArgumentParser: argument parsers.
    """
    # Define input and output field help message
    fields = dedent(
             '''
             output files:
                 lineages
                    successfully processed records.
                 lineages-fail
                    database records failed processing.

             required fields:
                 SEQUENCE_ID, SEQUENCE_INPUT, SEQUENCE_IMGT,
                 GERMLINE_IMGT_D_MASK, V_CALL, J_CALL
              ''')

    # Parent parser
    parser_parent = getCommonArgParser(db_out=False, log=True, format=False)

    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[parser_parent],
                            formatter_class=CommonHelpFormatter, add_help=False)

    group = parser.add_argument_group('tree building arguments')
    group.add_argument('--collapse', action='store_true', dest='collapse',
                        help='''Collapse identical sequences.''')
    group.add_argument('--md', nargs='+', action='store', dest='meta_data',default=None,
                       help='''Include metadata in sequence ID.''')

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main
    """
    # Parse command line arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    del args_dict['db_files']

    # Call main for each input file
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        buildTrees(**args_dict)
