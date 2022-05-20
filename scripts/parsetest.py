from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

#new 683.4873559474945
#old 621.9745631217957
#new 382.72539806365967
def parsefasta_new3(filename):
    with open(filename) as fh:
        sequence=""
        sequences=[]
        count=0
        for line in fh:
            if line[0] != '>':           
                sequence+=line.rstrip("\n")
                count+=1
            elif line[0] == '>' and count==0:
                continue
            elif line[0] == '>' and count>0:
                sequences.append(sequence)
                sequence=""
                count=0
        sequences.append(sequence)
        xmlfilename = "scriptxml_new3.xml"
        with open(xmlfilename, "w") as save_file:
            for seq in sequences:
                print(seq)    
                #result_handle = NCBIWWW.qblast("blastn", "nt", sequence=seq, format_type="XML") 
                #save_file.write(result_handle.read())
        return xmlfilename

def parsefasta_new(filename):
    xmlfilename = "scriptxml_new.xml"
    with open(filename) as fh:
        with open(xmlfilename, "w") as save_file:
             for seq_record in SimpleFastaParser(fh):
                print(seq_record[1])
                result_handle = NCBIWWW.qblast("blastn", "nt", sequence=seq_record[1], format_type="XML") 
                save_file.write(result_handle.read())
    return xmlfilename

def parsefasta_old(filename):
    xmlfilename = "scriptxml_old.xml"
    sequences = SeqIO.parse(open(filename,mode="r"),"fasta")
    with open(xmlfilename, "w") as save_file:
        for seq_record in sequences:
            print(seq_record.seq)
            #result_handle = NCBIWWW.qblast("blastn", "nt", sequence=seq_record.seq, format_type="XML") 
            #save_file.write(result_handle.read())
    return xmlfilename

def parsexml_old(filename):
    items=[]
    for record in NCBIXML.parse(open(filename)):
        if record.alignments:
            for align in record.alignments:
                for hsp in align.hsps: 
                    item=dict({"sequence":align.title,"length":align.length,"query":hsp.query[:50]})
                    items.append(item)
    return items

def parsexml(filename):
    items=[]
    recs=[]
    i = 0
    for record in NCBIXML.parse(open(filename)):
        i+=1
        if record.alignments:
            print(i)
            print("record.query_id: %s " % record.query_id)
            rec=dict({"query_id":record.query_id})
            recs.append(rec)
            for align in record.alignments:
                for hsp in align.hsps: 
                    item=dict({"iteration":i,"query_id":record.query_id,"query_length":record.query_length,"sequence":align.title,"length":align.length })
                    items.append(item)
    return items


##main
import time
start=time.time()
filename = "blast.fasta"
#xmlfile=parsefasta_new3(filename)
#xmlfile=parsefasta_old(filename)
xmlfile=parsefasta_new(filename)
#result=parsexml("scriptxml_new.xml")
#print(result)
end=time.time()
totaltime=end-start
print("\n"+str(totaltime))

