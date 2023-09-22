import os
import pandas as pd
import sys, getopt
import subprocess
import pprint

from collections import Counter
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

#dirname = "D:/PROJECT/Dongeubogam/tbc2/6.otu"
#ucfile = "D:/PROJECT/Dongeubogam/tbc2/5.chimera/all.tbc.nonchi.derep.uc"
#fastafile = "D:/PROJECT/Dongeubogam/tbc2/5.chimera/all.tbc.nonchi.derep.fasta"
#taxfile = "D:/PROJECT/Dongeubogam/tbc2/6.otu/ez/all.tbc.nonchi.derep_tax_assignments_w_tax.txt"


def subprocess_open(command, shell_use=True):
    popen = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell_use)
    (stdoutdata, stderrdata) = popen.communicate()    
    return stdoutdata, stderrdata

def executeVsearch(file, outdir, threads=4):
    if(outdir is None):
        outdir = os.path.dirname(file)
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    ucout = os.path.join(outdir,"un_cluster.uc")
    uctab = os.path.join(outdir,"un_otu.tab")
    command = "vsearch --threads "+str(threads)+" \
            --cluster_size "+file+" \
            --id 0.97 \
            --strand plus \
            --sizein \
            --sizeout \
            --uc "+ucout+" \
            --otutabout "+uctab
    #print(command)
    subprocess_open(command)
    return ucout, uctab
                              
def parseSize(line):    
    info = line.split(";size=")
    seqId = info[0]
    sample = info[0][0:info[0].rfind('.')]
    size = int(info[1][:-1])
    return [seqId,sample,size]
            
def parseUC(filename):
    rep_dict = defaultdict(lambda: Counter())
    
    seq_label_index = 8
    rep_label_index = 9
    
    with open(filename,"r") as f:
        line_count = 0
        for line in f:            
            line_count += 1
            if( line_count % 300000 == 0):
                print(line_count)
            line = line.strip()            
            if line.startswith('#') or len(line) == 0:
                continue
            if(line.startswith('S')):
                tmp = line.split('\t')
                tmpId = tmp[seq_label_index]
                info = parseSize(tmpId)
                c = Counter() 
                c[tmpId] = info[2]
                rep_dict[info[0]] = c
            elif(line.startswith('H')):
                tmp = line.split('\t')
                hitId = tmp[seq_label_index]
                hitInfo = parseSize(hitId)
                repId = tmp[rep_label_index]
                repInfo = parseSize(repId)                
                rep_dict[repInfo[0]][hitId] = hitInfo[2]
    return (rep_dict)




def pickCentroid(filename,tax_to_name,rep_dict,outdir):
    if(outdir is None):
        outdir = os.path.dirname(file)
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    bname = os.path.basename(filename)    
    outfile = os.path.join(outdir,bname.replace(".fasta",".tbc.centroids.fasta"))
    fw = open(outfile,"w")
    name_otu = defaultdict()
    for key,value in rep_dict.items():
        otu_id = tax_to_name[key]
        name_otu[value[0]] = otu_id
   
    with open(filename,"r") as fr:
        for title, seq in SimpleFastaParser(fr):
            if(name_otu.get(title) is not None):
                fw.write(">"+name_otu[title]+"\n")
                fw.write(seq+"\n")
    fw.close()
    return outfile

def pickUnassigned(filename,un_dict,outdir):
    if(outdir is None):
        outdir = os.path.dirname(file)
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    bname = os.path.basename(filename)    
    outfile = os.path.join(outdir,bname.replace(".fasta",".tbc.unassigned.fasta"))
    fw = open(outfile,"w")
#     un_dict_key = set(un_dict.keys())
    line_count = 0
    with open(filename,"r") as fr:
        for title, seq in SimpleFastaParser(fr):           
            line_count += 1   
            if( line_count % 100000 == 0):
                print(line_count)
            if(un_dict.get(title) is not None):
                for key,value in un_dict[title].items():                           
                    fw.write(">"+key+"\n")
                    fw.write(seq+"\n")
#                 un_dict_key.remove(title)

    fw.close()
    return outfile

def parseTax(taxfile,rep_dict):
    unassigned={}
    tax_dict = {}
    rep_tax = {}
    with open(taxfile, 'r') as f:
        line_count = 0
        for line in f:            
            line_count += 1
            if( line_count % 100000 == 0):
                print(line_count)
            contents = line.split("\t")
            tax = contents[1]

            info = parseSize(contents[0])
            full_name = contents[0]
            # if(tax.startswith("Unassigned")):
            #     ### unassigned[full_name]=rep_dict[info[0]]
            #     tmp_un_dict = {}
            #     tmp_un_dict[info[1]] = info[2]
            #     unassigned[full_name]= tmp_un_dict
            #     continue
            
            tmp_rep_dict = {}
            for key,value in rep_dict[info[0]].items():
                # print(key+"\t"+str(value))
                tmpInfo = parseSize(key)
                tmp_rep_dict[tmpInfo[1]] = tmpInfo[2]
            if(tax_dict.get(tax) is None):
                tax_dict[tax] = tmp_rep_dict
                rep_tax[tax] = [full_name,info[2]]
            else:
                repInfo = rep_tax[tax]
                if(repInfo[1] < info[2]):
                    rep_tax[tax] = [full_name,info[2]]
                
                tmp_tax_dict = tax_dict[tax]
                tmp_dict = tmp_rep_dict
                for currentkey,value in tmp_dict.items():                         
                    tmp_size = tmp_dict[currentkey]
                    if (tmp_tax_dict.get(currentkey) is None):
                        tax_dict[tax][currentkey] = tmp_size                        
                    else:
                        tax_dict[tax][currentkey]+= tmp_size

    return tax_dict, unassigned, rep_tax

def mergeResult(uctab,tax_dict):
    
    df = pd.DataFrame().from_dict(tax_dict,orient='index')
    df = df.fillna(0).astype(int)
    df['taxonomy'] = df.index

    new_index = []
    tax_to_name = {}
    for i in range(len(df.index)):
        otu_id = "OTU_"+str(i)
        tax_to_name[df.index[i]] = otu_id
        new_index.append(otu_id)
    df.index = new_index    
    df.index.name = "#OTU ID"
    print(df)
       
    # undf = pd.read_csv(uctab,header=0,sep="\t",index_col=0)
    # undf
    # undf['taxonomy']="Unassgined"
    # # print(df['taxonomy'])

    # finaldf = pd.concat([df,undf],sort=True).fillna(0,downcast=True)
    # # print(pd.concat([df,undf],sort=True))
    # print(undf['taxonomy'])


    # new_index = []
    # tax_to_name = {}
    # for i in range(len(finaldf.index)):
    #     otu_id = "OTU_"+str(i)
    #     tax_to_name[finaldf.index[i]] = otu_id
    #     new_index.append(otu_id)
    # finaldf.index = new_index    
    # finaldf.index.name = "#OTU ID"

    # mapout = os.path.join(os.path.dirname(uctab),"otu_table_w.txt")
    # finaldf.to_csv(mapout,sep="\t")
    
    # map_mc2_out = os.path.join(os.path.dirname(uctab),"otu_table_mc2_w.txt")
    # finaldf[finaldf.sum(axis=1) > 1].to_csv(map_mc2_out,sep="\t")

#     finaldf.to_csv((map_mc2_out,sep="\t")
    return tax_to_name
    

def execute(taxfile, fastafile, ucfile, threads, outdir):
    print("Parse UC")
    cluster = parseUC(ucfile)
    print("Parse Tax")
    tax_dict, unassigned , rep_tax = parseTax(taxfile,cluster)
    print("Parse Unassigned")
    unassigned = pickUnassigned(fastafile,unassigned,outdir)
    print("Execute Vsearch")
    ucout, uctab = executeVsearch(unassigned,outdir,threads)
    print("Merge Result")
    tax_to_name = mergeResult(uctab,tax_dict)    
    print("Pick Centroid")
    centroid = pickCentroid(fastafile,tax_to_name,rep_tax, outdir)


#outdir = "D:/PROJECT/Dongeubogam/tbc2/6.otu/"
#uctab = os.path.join(outdir,"un_otu.tab")

# cluster = parseUC(ucfile)
# tax_dict, unassigned , rep_tax = parseTax(taxfile,cluster)
# tax_to_name = mergeResult(uctab,tax_dict)    
# # print(unassigned.values())

def simpleParseUC(filename):
    rep_dict = defaultdict(lambda: Counter())
    seq_label_index = 8
    rep_label_index = 9
    id_set = set()
    with open(filename,"r") as f:
        line_count = 0
        for line in f:            
            line_count += 1
            if( line_count % 300000 == 0):
                print(line_count)
            line = line.strip()            
            if line.startswith('#') or len(line) == 0:
                continue
            if(line.startswith('S')):
                tmp = line.split('\t')
                tmpId = tmp[seq_label_index]
                info = parseSize(tmpId)
                c = Counter() 
                c[info[1]] = info[2]
                rep_dict[info[0]] = c
                id_set.add(info[1])
            elif(line.startswith('H')):
                tmp = line.split('\t')
                hitId = tmp[seq_label_index]
                hitInfo = parseSize(hitId)
                repId = tmp[rep_label_index]
                repInfo = parseSize(repId)                
                rep_dict[repInfo[0]][hitInfo[1]] = hitInfo[2]
                id_set.add(hitInfo[1])
    return rep_dict, sorted(id_set)

def simpleParseTax(taxfile):
    tax_dict = {}
    with open(taxfile, 'r') as f:
        line_count = 0
        for line in f:            
            line_count += 1
            if( line_count % 100000 == 0):
                print(line_count)
            contents = line.split("\t")
            info = parseSize(contents[0])
            tax_dict[info[0]] = contents[1]

    return tax_dict

def simpleASV(taxfile, fastafile, ucfile, outdir):
    tax_dict = simpleParseTax(taxfile)
    rep_dict, id_set = simpleParseUC(ucfile)
    
    id_set = sorted(id_set)
    id_dic = dict.fromkeys(id_set, 0) 
    for i,key in enumerate(id_dic.keys()):
        id_dic[key] = i
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    table_file = os.path.join(outdir,"otu_table_w.txt")    
    table_mc2_file = os.path.join(outdir,"otu_table_mc2_w.txt")
    
    with open(table_file,"w") as fw:
        header_string = "#OTU ID\t"+"\t".join(str(x) for x in id_dic.keys())+"\ttaxonomy\n"
        fw.write(header_string)
        with open(table_mc2_file,"w") as fw2:
            fw2.write(header_string)
            for key,inner_dict in rep_dict.items():
                tmparray = [0] * len(id_dic)
                for inner_key in inner_dict.keys():
                    arridx = id_dic[inner_key]
                    tmparray[arridx] = inner_dict[inner_key]
                one_line = key+"\t"+"\t".join(str(x) for x in tmparray)+"\t"+tax_dict[key]+"\n"
                fw.write(one_line)
                if sum(tmparray) > 1:
                    fw2.write(one_line)

#simpleASV(taxfile, fastafile, ucfile, outdir)

def usage():
    print("make_otu_table.py -t <taxfile> -u <ucfile> -f <fastafile> -c <threads> -o <outdir> -m <method>")
                    
def main():
    taxfile=""
    ucfile=""
    fastafile=""
    outdir = None
    ncore = 1
    method = 1
    try:
        opts, args = getopt.getopt(sys.argv[1:],"ht:u:f:c:o:m:")
    except getopt.GetoptError as err: 
        usage()
        sys.exit(1)
    
    if len(sys.argv)<2:
        usage()
        return
    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            return
        elif opt in ("-t"):
            taxfile = arg
        elif opt in ("-u"):
            ucfile = arg
        elif opt in ("-c"):
            ncore = arg
        elif opt in ("-f"):
            fastafile = arg
        elif opt in ("-o"):
            outdir = arg
        elif opt in ("-m"):
            method = arg
        else:
            usage()
    simpleASV(taxfile, fastafile, ucfile, outdir)
    if method == 2:
        print(method)      
#            execute(taxfile, fastafile, ucfile, ncore, outdir)
       
if __name__ == "__main__": 
    main()




