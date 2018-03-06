#!/usr/bin/env python3

import changeo.Parsers

def outputIgPhyML(clones, sequences, outdir):
    """
    1. Remove columns in clonal lineage that only contain gaps in data seqs (ignore reference)
    2. Remove same columns in refernce
    3. Create list of IMGT numbering for each site in alignment
    :param
    clones: distrionary of clones
    sequences: list of masked sequences
    outdir: string specifying output directory
    :return: void, outputs files to "outdir"
    """
    sites = len(sequences[0])
    for i in sequences:
        if len(i) != sites:
            print("Sequences within clone %d are not the same length!" % clones[0].clone)
            exit(1)

    nseqs = len(sequences)
    germline = clones[0].getField("germline_imgt_d_mask")
    if len(sequences[0]) % 3 != 0:
        print("number of sites must be divisible by 3! %d", len(sequences))
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

    for i in range(0, sites, 3):
        if tallies[i//3] > 0:
            newgerm.append(germline[i:(i+3)])
            imgt.append(i//3)

    # Output fasta file of masked, concatenated sequences
    outfile = outdir+"/"+clones[0].clone+".fa"
    clonef = open(outfile, 'w')
    for j in range(0,nseqs):
        print(">%s" % clones[j].sequence_id,file=clonef)
        for i in range(0,len(newgerm)):
            print("%s" % newseqs[j][i],end='',file=clonef)
        print("\n",end='',file=clonef)
    print(">%s_GERM" % clones[0].clone, file=clonef)
    for i in range(0, len(newgerm)):
        print("%s" % newgerm[i], end='', file=clonef)
    print("\n", end='', file=clonef)
    clonef.close()

    #output partition file
    partfile = outdir + "/" + clones[0].clone + ".part.txt"
    partf = open(partfile, 'w')
    print("%d %d" % (2, len(newgerm)), file=partf)
    print("FWR:IMGT", file=partf)
    print("CDR:IMGT", file=partf)
    print("%s" % (clones[0].v_call.split("*")[0]),file=partf)
    print("%s" % (clones[0].j_call.split("*")[0]),file=partf)
    print(",".join(map(str,imgt)), file=partf)


handle = open('test_proc/test_db-pass_clone-pass_germ-pass.tsv')
records = changeo.Parsers.ChangeoReader(handle)

cloneseqs = {}
cloneids = {}
clonegerms = {}
clones = {}
outdir = "test_proc/igphyml_out"
for r in records:
    if r.functional == True:
        mask_seq = changeo.Parsers.maskSplitCodons(r)
        if r.clone in clones:
            clones[r.clone].append(r)
            cloneseqs[r.clone].append(mask_seq)
        else:
            clones[r.clone]=[r]
            cloneseqs[r.clone] = [mask_seq]
    else:
        print("Skipping %s", r.sequence_id)

clonesizes={}
for k in clones.keys():
    changeo.Parsers.outputIgPhyML(clones[str(k)], cloneseqs[str(k)], outdir)
    clonesizes[str(k)] = len(cloneseqs[str(k)])

repfile1 = outdir + "/repFile.tsv"
repfile2 = outdir + "/repFile.N.tsv"
rout1 = open(repfile1, 'w')
rout2 = open(repfile2, 'w')
print(len(clonesizes),file=rout1)
print(len(clonesizes),file=rout2)
for key in sorted(clonesizes, key=clonesizes.get,reverse=True):
    print(key + "\t"+str(clonesizes[key]))
    outfile = outdir + "/" + key + ".fa"
    partfile = outdir + "/" + key + ".part.txt"
    print("%s\t%s\t%s\t%s" % (outfile, "N", key+"_GERM", partfile), file=rout1)
    print("%s\t%s\t%s\t%s" % (outfile, "N", key+"_GERM", "N"), file=rout2)

handle.close()