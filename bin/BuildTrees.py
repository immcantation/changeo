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
from changeo.IO import AIRRReader, ChangeoReader, AIRRWriter, ChangeoWriter, splitFileName, getDbFields
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

    # adjust starting position of query sequence
    qi = qi[(receptor.v_seq_start - 1):]

    # deal with the fact that it's possible to start mid-codon
    scodons = [si[i:i + 3] for i in range(0, len(si), 3)]
    for i in range(0, len(scodons)):
        if scodons[i] != '...':
            if scodons[i][0:2] == '..':
                scodons[i] = "NN"+scodons[i][2]
                qi = "NN" + qi
                spos = i
                break
            elif scodons[i][0] == '.':
                scodons[i] = "N" + scodons[i][1:3]
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
                #print("But codon was apparently preserved")
                #log[str(spos)]="IN-FRAME"
                if 'IN-FRAME' in  log:
                    log['IN-FRAME'] = log['IN-FRAME'] + "," +  str(spos)
                else:
                    log['IN-FRAME'] = str(spos)
            elif qpos >= len(qcodons):
                log['PASS'] = False
                log['FAIL'] = "FAILED_MATCH_QSTRING:"+str(spos)
            elif spos >= len(scodons) or qcodons[qpos] != scodons[spos]:
                scodons[ospos] = "NNN"
                if spos >= len(scodons):
                    #print("Masked %s at position %d, at end of subject sequence" % (scodons[ospos], ospos))
                    #log[str(spos)] = "END"
                    if 'END-MASKED' in log:
                        log['END-MASKED'] = log['END-MASKED'] + "," + str(spos)
                    else:
                        log['END-MASKED'] = str(spos)
                else:
                    #print("Masked %s at position %d, but couldn't find upstream match" % (scodons[ospos], ospos))
                    #log[str(spos)] = "FAILED_MATCH"
                    log['PASS']=False
                    log['FAIL']="FAILED_MATCH:"+str(spos)
                    #exit(1)
            elif qcodons[qpos] == scodons[spos]:
                #print("Masked %s at position %d" % (scodons[ospos], ospos))
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

def outputIgPhyML(clones, sequences, collapse=False, logs=None, fail_writer=None, out_dir=None):
    """
    Create intermediate sequence alignment and partition files for IgPhyML output

    Arguments:
        clones (list): receptor objects within the same clone.
        sequences (list): sequences within the same clone (share indexes with clones parameter).
        collapse (bool): if True collapse identical sequences.
        logs (dict of OrderedDicts): contains log information for each sequence
        out_dir (str): directory for output files.
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


    useqs = {}
    conseqs = []
    for j in range(0, nseqs):
        conseq = "".join([str(seq_rec) for seq_rec in newseqs[j]])
        if conseq in useqs and collapse:
            logs[clones[j].sequence_id]['PASS'] = False
            logs[clones[j].sequence_id]['FAIL'] = "Duplication of " + clones[useqs[conseq]].sequence_id
            logs[clones[j].sequence_id]['DUPLICATE']=True
            if fail_writer is not None:
                fail_writer.writeReceptor(clones[j])
        else:
            useqs[conseq] = j
        conseqs.append(conseq)

    # Output fasta file of masked, concatenated sequences
    outfile = out_dir + "/" + clones[0].clone + ".fa"
    clonef = open(outfile, 'w')
    if collapse:
        for seq, num in useqs.items():
            sid = clones[num].sequence_id.translate(transtable)
            print(">%s\n%s" % (sid, seq), file=clonef)
            if len(useqs) == 1 and duplicate:
                print(">%s_1\n%s" % (sid, seq), file=clonef)
    else:
        for j in range(0, nseqs):
            sid = clones[j].sequence_id.translate(transtable)
            print(">%s\n%s" % (sid, conseqs[j]), file=clonef)
        if nseqs == 1 and duplicate:
            sid = clones[j].sequence_id.translate(transtable)
            print(">%s_1\n%s" % (sid, conseqs[j]), file=clonef)

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
        return len(useqs)
    else:
        return nseqs


def buildTrees(db_file, collapse=False, format=default_format, out_args=default_out_args):
    """
    Masks codons split by alignment to IMGT reference, then produces input files for IgPhyML

    Arguments:
        db_file (str): input tab-delimited database file.
        collapse (bool): if True collapse identical sequences.
        format (str): input and output format.
        out_args (dict): arguments for output preferences.
    Returns:
        None: outputs files for IgPhyML.

    """
    start_time = time()
    # TODO: changeo console log
    # pass_handle.name = "outdir/out_name_lineages.tsv"
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
    if format == 'changeo':
        reader = ChangeoReader
        writer = ChangeoWriter
        out_fields = getDbFields(db_file, reader=reader)
    elif format == 'airr':
        reader = AIRRReader
        writer = AIRRWriter
        out_fields = getDbFields(db_file, reader=reader)
    else:
        sys.exit('Error:  Invalid format %s' % format)

    # open input file
    handle = open(db_file, 'r')
    records = reader(handle)

    # TODO: make consistent with changeo behavior regarding --log and --failed
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
            log['FAIL'] = "NONFUNCTIONAL"
            log['SEQ_IN'] = r.sequence_input
            logs[r.sequence_id]=log
            if out_args['failed']:
                fail_writer.writeReceptor(r)
            seq_fail += 1

    clonesizes = {}
    pass_count = 0
    for k in clones.keys():
        clonesizes[str(k)] = outputIgPhyML(clones[str(k)], cloneseqs[str(k)], collapse=collapse,
                                           logs=logs,fail_writer=fail_writer,out_dir=clone_dir)
        pass_count += clonesizes[str(k)]
    fail_count = rec_count - pass_count

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
        # TODO: use file.write()
        print("%s\t%s\t%s\t%s" % (outfile, "N", key+"_GERM", partfile), file=pass_handle)

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
                 germ-pass
                    database with assigned germline sequences.
                 germ-fail
                    database with records failing germline assignment.

             required fields:
                 SEQUENCE_ID, SEQUENCE_VDJ or SEQUENCE_IMGT,
                 V_CALL or V_CALL_GENOTYPED, D_CALL, J_CALL,
                 V_SEQ_START, V_SEQ_LENGTH, V_GERM_START_IMGT, V_GERM_LENGTH_IMGT,
                 D_SEQ_START, D_SEQ_LENGTH, D_GERM_START, D_GERM_LENGTH,
                 J_SEQ_START, J_SEQ_LENGTH, J_GERM_START, J_GERM_LENGTH,
                 NP1_LENGTH, NP2_LENGTH

             optional fields:
                 N1_LENGTH, N2_LENGTH, P3V_LENGTH, P5D_LENGTH, P3D_LENGTH, P5J_LENGTH,
                 CLONE


             output fields:
                 GERMLINE_VDJ, GERMLINE_VDJ_D_MASK, GERMLINE_VDJ_V_REGION,
                 GERMLINE_IMGT, GERMLINE_IMGT_D_MASK, GERMLINE_IMGT_V_REGION,
                 GERMLINE_V_CALL, GERMLINE_D_CALL, GERMLINE_J_CALL,
                 GERMLINE_REGIONS
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
