from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys, getopt


def bokulich(infile,outfile, min_len=350, max_len=550):
    r = 3
    p = 0.75
    n = 0
    q = 3
#    print(min_len, max_len)
    
    with open(infile,"r") as in_handle:
        with open(outfile,"w") as fw:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                L = len(seq)
                BadRun = 0
                Ns = 0
                M = L
                for i in range(0, L):
                    Q = ord(qual[i]) - 33
                    if Q > q:
                        BadRun = 0
                    else:
                        BadRun += 1
                        if BadRun > r:
                            M = i - BadRun + 1
                            if float(M)/L < p:
                                M = 0
                            break

                    if seq[i] == 'N':
                        Ns += 1
                        if Ns > n:
                            M = 0
                            break
                new_len = len(seq[0:M]) 
                if( M > 0 and new_len >= min_len and new_len <= max_len):
                    fw.write("@"+title+"\n")
                    fw.write(seq[0:M]+"\n")
                    fw.write("+"+"\n")
                    fw.write(qual[0:M]+"\n")
             
                    
def usage():
    print("bokulich.py -i <infilename> -o <outfilename> -m <min_len> -M <max_len>")
                    
def main():
    infilename=""
    outfilename=""
    min_len = 350
    max_len = 550
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:m:M:",["in=","out="])
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
        elif opt in ("-i", "--in"):
            infilename = arg
        elif opt in ("-o", "--out"):
            outfilename = arg
        elif opt in ("-m"):
            min_len = int(arg)
        elif opt in ("-M"):
            max_len = int(arg)
        else:
            usage()
        
    bokulich(infilename,outfilename,min_len,max_len)
       
if __name__ == "__main__": 
    main()

        

#bokulich("D:/PROJECT/GS001_1.fastq","D:/PROJECT/GS001_2.fastq",10,500)

