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

def maskSplitCodons(receptor,recursive=False):
    """
    Identify junction region by IMGT definition.

    Arguments:
      receptor : Receptor object.
      recursive (bool) : was this method part of a recursive call?
    Returns:
      str : modified IMGT gapped sequence.
    """
    debug = False
    qi = receptor.sequence_input
    si = receptor.sequence_imgt
    log = OrderedDict()
    log['ID']=receptor.sequence_id
    log['CLONE']=receptor.clone
    log['PASS'] = True
    if debug:
        print(receptor.sequence_id)
    # adjust starting position of query sequence
    qi = qi[(receptor.v_seq_start - 1):]

    #TODO: tally where --- gaps are in IMGT sequence and remove them for now
    gaps = []
    nsi = ""
    for i in range(0,len(si)):
        if si[i] == "-":
            gaps.append(1)
        else:
            gaps.append(0)
            nsi = nsi + si[i]

    #TODO: find any gaps not divisble by three
    curgap = 0
    for i in gaps:
        if i == 1:
            curgap += 1
        elif i == 0 and curgap != 0:
            if curgap % 3 != 0 :
                if debug:
                    print("Frame-shifting gap detected! Cowardly refusing to include sequence.")
                log['PASS'] = False
                log['FAIL'] = "FRAME-SHIFTING DELETION"
                log['SEQ_IN'] = receptor.sequence_input
                log['SEQ_IMGT'] = receptor.sequence_imgt
                log["SEQ_MASKED"] = receptor.sequence_imgt
                return receptor.sequence_imgt, log
            else:
                curgap = 0

    si = nsi

    # deal with the fact that it's possible to start mid-codon
    scodons = [si[i:i + 3] for i in range(0, len(si), 3)]
    for i in range(0, len(scodons)):
        if debug:
            print(scodons[i],qi[0:3])
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

    frameshifts = 0
    s_end = 0 #adjust for the fact that IMGT sequences can end on gaps
    for i in range(spos, len(scodons)):
        if scodons[i] != '...' and len(scodons[i]) == 3:
            s_end = i
    if debug:
        print(str(s_end) + ":" + str(len(scodons)))
        print(scodons[s_end])
    s_end += 1
    qpos = 0
    # TODO: for loop with zip()
    while spos < s_end and qpos < len(qcodons):
        if debug:
            print(scodons[spos] + "\t" + qcodons[qpos])
        if scodons[spos] == '...' and qcodons[qpos] != '...': #if IMGT gap, move forward in imgt
            spos += 1
        elif scodons[spos] == qcodons[qpos]: # if both are the same, move both forward
            spos += 1
            qpos += 1
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
                if debug:
                    print("Checking " + scodons[spos]+ "\t" + qcodons[qpos])
                qpos += 1
            if qcodons[qpos-1] == scodons[ospos]: #if codon in previous position is equal to original codon, it was preserved
                qpos -= 1
                spos = ospos
                if debug:
                    print("But codon was apparently preserved")
                if 'IN-FRAME' in  log:
                    log['IN-FRAME'] = log['IN-FRAME'] + "," +  str(spos)
                else:
                    log['IN-FRAME'] = str(spos)
            elif qpos >= len(qcodons) and spos < s_end:
                if debug:
                    print("FAILING MATCH")
                log['PASS'] = False #if no match for the adjacent codon was found, something's up.
                log['FAIL'] = "FAILED_MATCH_QSTRING:"+str(spos)
                #figure out if this was due to a frame-shift by repeating this method but with an edited input sequence
                if not recursive:
                    for ins in range(1,3):
                        ros = receptor.sequence_input
                        ris = receptor.sequence_imgt
                        psite = receptor.v_seq_start-1+oqpos*3
                        pisite = ospos * 3
                        if (psite+3 + ins) < len(ros) and (pisite+3) < len(ris):
                        #cut out 1 or 2 nucleotides downstream of offending codon
                            receptor.sequence_input = ros[0:(psite+3)]+ros[(psite+3 + ins):]
                            receptor.sequence_imgt = ris[0:(pisite+3)] + ris[(pisite+3):]
                            if debug:
                                print(ros + "\n"+receptor.sequence_input)
                                print(ris + "\n" + receptor.sequence_imgt)
                                print("RUNNING %d\n"%ins)
                            mout = maskSplitCodons(receptor,recursive=True)
                            if mout[1]['PASS']:
                                #if debug:
                                receptor.sequence_input = ros
                                receptor.sequence_imgt = ris
                                frameshifts += 1
                                if debug:
                                    print("FRAMESHIFT of length %d!" % ins)
                                log['FAIL'] = "SINGLE FRAME-SHIFTING INSERTION"
                                break
                            else:
                                receptor.sequence_input = ros
                                receptor.sequence_imgt = ris

                            '''
                                message,site = mout[1]['FAIL'].split(":")
                                if message == "FAILED_MATCH_QSTRING":
                                    log['FAIL'] = "FRAME-SHIFTING INSERTION"
                                    site = int(site)
                                    if site != spos:
                                        if debug:
                                            print("frame fixed site %d %d" % (spos, site))
                                        frameshifts += 1
                                        qi = receptor.sequence_input
                                        qi = qi[(receptor.v_seq_start - 1):]
                                        si = receptor.sequence_imgt
                                        nsi = ""
                                        for i in range(0, len(si)):
                                            if si[i] != "-":
                                               nsi = nsi + si[i]
                                        si = nsi
                                        scodons = [si[i:i + 3] for i in range(0, len(si), 3)]
                                        qcodons = [qi[i:i + 3] for i in range(0, len(qi), 3)]
                                        spos = ospos
                                        qpos = oqpos
                                        if debug:
                                            print("just re-assigned "+scodons[spos] + "\t" + qcodons[qpos])
                                        break
                                if debug:
                                    print(message)
                                receptor.sequence_input = ros
                                receptor.sequence_imgt = ris
                                log['FAIL'] = "FAILED_MATCH_QSTRING:"+str(spos)
                                    #print(mout[1]['FAIL'] + "\t" + str(qpos) + "\t" + str(spos))
                            '''
            elif spos >= s_end or qcodons[qpos] != scodons[spos]:
                scodons[ospos] = "NNN"
                if spos >= s_end:
                    if debug:
                        print("Masked %s at position %d, at end of subject sequence" % (scodons[ospos], ospos))
                    #log[str(spos)] = "END"
                    if 'END-MASKED' in log:
                        log['END-MASKED'] = log['END-MASKED'] + "," + str(spos)
                    else:
                        log['END-MASKED'] = str(spos)
                else:
                    if debug:
                        print("Masked %s at position %d, but couldn't find upstream match" % (scodons[ospos], ospos))
                    log['PASS']=False
                    log['FAIL']="FAILED_MATCH:"+str(spos)
            elif qcodons[qpos] == scodons[spos]:
                if debug:
                    print("Masked %s at position %d" % (scodons[ospos], ospos))
                scodons[ospos] = "NNN"
                if 'MASKED' in  log:
                    log['MASKED'] = log['MASKED'] + "," + str(spos)
                else:
                    log['MASKED'] = str(spos)

            else:
                log['PASS'] = False
                log['FAIL'] = "UNKNOWN"

    if not log['PASS'] and not recursive:
        #if log['FAIL'] == "FRAME-SHIFTING INSERTION":
        log['FRAMESHIFTS'] = frameshifts
    if len(scodons[-1]) != 3:
        if scodons[-1] == ".." or scodons[-1] == ".":
            scodons[-1] = "..."
        else:
            scodons[-1] = "NNN"
        #log[str(len(scodons))] = "MASKED"
        if 'MASKED' in log:
            log['MASKED'] = log['MASKED'] + "," + str(len(scodons))
        else:
            log['MASKED'] = str(spos)

    concatenated_seq = Seq("")
    for i in scodons:
        concatenated_seq += i

    # TODO: add --- gaps back to IMGT sequence
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
    log['SEQ_IN'] = receptor.sequence_input
    log['SEQ_IMGT'] = receptor.sequence_imgt
    log["SEQ_MASKED"] = concatenated_seq

    return concatenated_seq, log

def outputIgPhyML(clones, sequences, meta_data=None, collapse=False, logs=None, fail_writer=None, out_dir=None,
                  min_seq=0):
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
        min_seq (int): minimum number of data sequences to include
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
        print("number of sites must be divisible by 3! len: %d, clone: %s , seq: %s" %(len(sequences[0]),clones[0].clone,sequences[0]))
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
    #print("MINSEQ "+str(min_seq))
    if collapse and len(useqs_f) < min_seq:
        for seq_f, num in useqs_f.items():
            logs[clones[num].sequence_id]['FAIL'] = "Clone too small: "+ str(len(useqs_f))
            logs[clones[num].sequence_id]['PASS'] = False
        return -len(useqs_f);

    if not collapse and len(conseqs) < min_seq:
        for j in range(0, nseqs):
            logs[clones[j].sequence_id]['FAIL'] = "Clone too small: "+str(len(conseqs))
            logs[clones[j].sequence_id]['PASS'] = False
        return -len(conseqs);

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


def buildTrees(db_file, meta_data=None, collapse=False, min_seq=None, format=default_format, out_args=default_out_args):
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

    if min_seq == 0:
        min_seq = 0
    else:
        min_seq = int(min_seq[0])

    cloneseqs = {}
    clones = {}
    logs = OrderedDict()
    rec_count = 0
    seq_fail = 0
    nf_fail = 0
    del_fail = 0
    in_fail = 0
    minseq_fail = 0
    other_fail = 0
    for r in records:
        rec_count += 1
        #printProgress(rec_count, rec_count, 0.05, start_time)
        if r.functional:
            mout = maskSplitCodons(r)
            #Sometimes v start site is off a little
            '''if not mout[1]['PASS'] and r.v_seq_start > 0:
                r.v_seq_start -= 3
                mout2 = maskSplitCodons(r)
                r.v_seq_start += 3
                if mout2[1]['PASS']:
                    mout = mout2'''
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
                if mout[1]['FAIL'] == "FRAME-SHIFTING DELETION":
                    del_fail += 1
                elif mout[1]['FAIL'] == "SINGLE FRAME-SHIFTING INSERTION":
                    in_fail += 1
                else:
                    other_fail += 1

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
            nf_fail += 1

    clonesizes = {}
    pass_count = 0
    nclones = 0
    for k in clones.keys():
        clonesizes[str(k)] = outputIgPhyML(clones[str(k)], cloneseqs[str(k)], meta_data=meta_data, collapse=collapse,
                                           logs=logs,fail_writer=fail_writer,out_dir=clone_dir,min_seq=min_seq)
        #If clone is too small, size is returned as a negative
        if clonesizes[str(k)] > 0:
            nclones += 1
            pass_count += clonesizes[str(k)]
        else:
            seq_fail -= clonesizes[str(k)]
            minseq_fail  -= clonesizes[str(k)]
    fail_count = rec_count - pass_count

    log_handle = None
    if out_args['log_file'] is not None:
        log_handle = open(out_args['log_file'], 'w')
        for j in logs.keys():
            printLog(logs[j], handle=log_handle)

    # TODO: changeo console log
    print(nclones, file=pass_handle)
    for key in sorted(clonesizes, key=clonesizes.get, reverse=True):
        #print(key + "\t" + str(clonesizes[key]))
        outfile = clone_dir + "/" + key + ".fa"
        partfile = clone_dir + "/" + key + ".part.txt"
        if clonesizes[key] > 0:
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
    #log['TOTAL MASK FAIL'] = seq_fail
    log['NONFUNCTIONAL'] = nf_fail
    log['FRAMESHIFT_DEL'] = del_fail
    log['FRAMESHIFT_INS'] = in_fail
    log['CLONETOOSMALL'] = minseq_fail
    log['OTHER_FAIL'] = other_fail

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

    group.add_argument('--minseq', nargs='+', action='store', dest='min_seq', default=0,
                       help='''Minimum number of data sequences.''')

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
