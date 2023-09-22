import sys, getopt
import os.path

def tax2OTU(otutaxfile,reftaxfile,outfile):
    dict_tax = {}
    with open(reftaxfile, 'r') as f:
        for line in f:
            contents = line.strip().split("\t")
            sample =  contents[0].strip()
            tax = contents[1].strip()
            dict_tax[sample]=tax
    with open(otutaxfile, 'r') as f:
        with open(outfile, 'w') as fw:
            for line in f:
                contents = line.strip().split()
                if len(contents) == 4:
                    sample =  contents[1].strip()
                    if(sample in dict_tax):
                        tax =  dict_tax[sample]
                    else:
                        tax = sample
                    line = contents[0]+"\t"+tax+"\t"+contents[2]+"\t"+contents[3]
                fw.write(line+"\n")

def usage():
    print("id_to_fulltax.py -t <otutaxfile> -r <reftaxfile> -o <outfile> ")
                    
def main():
    otutaxfile=""
    reftaxfile=""
    outfile=""
    min_count=2
    try:
        opts, args = getopt.getopt(sys.argv[1:],"ht:r:o:",["in=","out="])
    except getopt.GetoptError as err: 
        usage()
        sys.exit(1)
    
    if len(sys.argv) < 2:
        usage()
        return
    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            return
        elif opt in ("-t"):
            otutaxfile = arg
        elif opt in ("-r"):
            reftaxfile = arg
        elif opt in ("-o"):
            outfile = arg
        else:
            usage()
            
    tax2OTU(otutaxfile,reftaxfile,outfile)
       
if __name__ == "__main__": 
    main()


