#!/usr/bin/env python3

from argparse import ArgumentParser
from textwrap import dedent

import sys
from time import time

# Presto and change imports
from presto.Defaults import default_out_args
from presto.IO import  printLog
from changeo.Defaults import default_format

import changeo
import changeo.Parsers
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs


def outputIgPhyML(clones, sequences, out_dir):
    """
    Create intermediate sequence alignment and partition files for IgPhyML output

    Arguments:
        clones (list): receptor objects within the same clone
        sequences (list): sequences within the same clone (share indexes with clones parameter)
        out_dir (string): directory for output files
    Returns:
        None. Outputs alignment and partition files
    """
    t=False
    duplicate = True # duplicate sequences in clones with only 1 sequence?
    sites = len(sequences[0])
    s=""
   # transtable = s.maketrans('.', '-')
    for i in sequences:
        if len(i) != sites:
            print("Sequences within clone %s are not the same length!" % clones[0].clone)
            exit(1)
       # i.translate(transtable)

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
        # print(newgerm[-1])
    elif len(lcodon) == 1:
        newgerm[-1] = newgerm[-1] + "NN"
        # print(newgerm[-1])


    # Output fasta file of masked, concatenated sequences
    outfile = out_dir+"/"+clones[0].clone+".fa"
    clonef = open(outfile, 'w')
    for j in range(0,nseqs):
        print(">%s" % clones[j].sequence_id,file=clonef)
        for i in range(0,len(newgerm)):
            print("%s" % newseqs[j][i],end='',file=clonef)
        print("\n",end='',file=clonef)
    if nseqs == 1 and duplicate:
        print(">%s" % clones[j].sequence_id+"_1", file=clonef)
        for i in range(0, len(newgerm)):
            print("%s" % newseqs[j][i], end='', file=clonef)
        print("\n", end='', file=clonef)
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
    print(",".join(map(str,imgt)), file=partf)



def buildTrees(db_file, out_args=default_out_args,format=default_format):
    """
    Masks codons split by alignment to IMGT reference, then produces input files for IgPhyML

    Arguments:
        db_file : input tab-delimited database file
        out_dir (string): directory to contain output files
        out_args: unused
        format: ununsed
    Returns:
        Nothing. Outputs files for IgPhyML

    """

    # Format options
    if format == 'changeo':
        reader = changeo.Parsers.ChangeoReader
    elif format == 'airr':
        reader = changeo.Parsers.AIRRReader
    else:
        sys.exit('Error:  Invalid format %s' % format)

    handle = open(db_file,'r')
    records = reader(handle) # open input file
    out_dir = out_args['out_dir'] # get output directory
    log_handle = open(out_args['log_file'],'w')

    failed = open(out_dir+"/"+"failedSeqs.fa",'w')

    cloneseqs = {}
    clones = {}
    for r in records:
        if r.functional:
            mout = changeo.Parsers.maskSplitCodons(r)
            mask_seq=mout[0]

            if mout[1]['PASS']:
                printLog(mout[1], handle=log_handle)
                if r.clone in clones:
                    clones[r.clone].append(r)
                    cloneseqs[r.clone].append(mask_seq)
                else:
                    clones[r.clone]=[r]
                    cloneseqs[r.clone] = [mask_seq]
            else:
                printLog(mout[1],handle=failed)
        else:
            print("Skipping %s", r.sequence_id)

    clonesizes={}
    for k in clones.keys():
        outputIgPhyML(clones[str(k)], cloneseqs[str(k)], out_dir)
        clonesizes[str(k)] = len(cloneseqs[str(k)])

    repfile1 = out_dir + "/repFile.tsv"
    repfile2 = out_dir + "/repFile.N.tsv"
    rout1 = open(repfile1, 'w')
    rout2 = open(repfile2, 'w')
    print(len(clonesizes),file=rout1)
    print(len(clonesizes),file=rout2)
    for key in sorted(clonesizes, key=clonesizes.get,reverse=True):
        print(key + "\t"+str(clonesizes[key]))
        outfile = out_dir + "/" + key + ".fa"
        partfile = out_dir + "/" + key + ".part.txt"
        print("%s\t%s\t%s\t%s" % (outfile, "N", key+"_GERM", partfile), file=rout1)
        print("%s\t%s\t%s\t%s" % (outfile, "N", key+"_GERM", "N"), file=rout2)

    handle.close()


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments:
    None

    Returns:
    an ArgumentParser object
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
    parser_parent = getCommonArgParser(log=True, format=True)

    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[parser_parent],
                            formatter_class=CommonHelpFormatter)

    # parser.add_argument('--version', action='version',
    #                    version='%(prog)s:' + ' %s-%s' %(__version__, __date__))

    # parser.add_argument('-r', nargs='+', action='store', dest='repo', required=True,
    #                   help='''List of folders and/or fasta files (with .fasta, .fna or .fa
    #                     extension) with germline sequences.''')

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
    #print(args_dict['out_dir'])

    # Call main for each input file
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        buildTrees(**args_dict)
