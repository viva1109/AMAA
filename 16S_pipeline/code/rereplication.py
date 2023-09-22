from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys, getopt
import re


def rereplicate(infile,outfile):

    with open(infile) as in_handle:        
        with open(outfile,"w") as fw:
            oldId = ""
            sample_count = 1;
            for title, seq in SimpleFastaParser(in_handle):
                tmp = title.split(";")
                sampleId = re.sub("\\.[0-9]*$","",tmp[0])
             #   print(sampleId)

                if(oldId != sampleId):
                    oldId = sampleId
                    sample_count = 1
                size = re.sub("size=","",tmp[1])
            #    print(size)
                for i in range(0,int(size)):
                    title = sampleId+"_"+str(sample_count)
                    sample_count+=1
                    fw.write(">"+title+"\n")
                    fw.write(seq+"\n")
                    
                
                    
def usage():
    print("rereplicate.py -i <infilename> -o <outfilename> ")
                    
def main():
    infilename=""
    outfilename=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["in=","out="])
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
        else:
            usage()
        
    rereplicate(infilename,outfilename)
       
if __name__ == "__main__": 
    main()

        

#bokulich("D:/PROJECT/GS001_1.fastq","D:/PROJECT/GS001_2.fastq",10,500)

#     rereplicate(infilename,outfilename)
