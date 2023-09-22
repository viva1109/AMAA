import sys, getopt
def tax2OTU(taxfile,otutabfile,outfile):
    dict_tax = {}
    with open(taxfile, 'r') as f:
        for line in f:
            contents = line.strip().split("\t")
            sample =  contents[0].strip()
            tax = contents[1].strip()
            dict_tax[sample]=tax
    with open(otutabfile, 'r') as f:
        with open(outfile, 'w') as fw:
            for line in f:
                if not line.startswith("#"):
                    contents = line.strip().split()
                    sample =  contents[0].strip()
                    tax =  dict_tax[sample]
                    line = line.strip()+"\t"+tax
                else:
                    line = line.strip()+"\ttaxonomy"                
                fw.write(line+"\n")
            
                    
def usage():
    print("match_otu_to_tax.py -t <taxfile> -m <otutafile> -o <outfile>")
                    
def main():
    taxfile=""
    otutafile=""
    outfile=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"ht:m:o:",["in=","out="])

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
        elif opt in ("-m"):
            otutafile = arg
        elif opt in ("-o"):
            outfile = arg
        else:
            usage()
        
    tax2OTU(taxfile,otutafile,outfile)
       
if __name__ == "__main__": 
    main()
