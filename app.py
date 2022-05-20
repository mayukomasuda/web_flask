from ast import parse
from fileinput import close
import os
from flask import Flask, render_template, request,redirect, url_for
from flask import send_from_directory
from werkzeug.utils import secure_filename
from werkzeug.datastructures import  FileStorage
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

app = Flask(__name__)   
app.config['UPLOAD_PATH'] = 'uploads'

def validate(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # Return True if this is fasta file format

def parsefasta(filename):
    xmlfilename = "results.xml"
    filepath = os.path.join("./",app.config['UPLOAD_PATH'], filename) 
    with open(filepath) as fh:
        with open(xmlfilename, "w") as save_file:
             for seq_record in SimpleFastaParser(fh):
                print(seq_record[1])
                result_handle = NCBIWWW.qblast("blastn", "nt", sequence=seq_record[1], format_type="XML") 
                save_file.write(result_handle.read())
    return xmlfilename

def parsexml(filename):
    items=[]
    i = 0
    for record in NCBIXML.parse(open(filename)):
        i+=1
        if record.alignments:
            for align in record.alignments:
                for hsp in align.hsps: 
                    item=dict({"iteration":i,"query_id":record.query_id,"query_length":record.query_length,"title":align.title,"length":align.length })
                    items.append(item)
    return items

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/', methods=['POST'])
def upload_files():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    filepath = os.path.join("./",app.config['UPLOAD_PATH'], filename) 
    if filename != '':
        uploaded_file.save(os.path.join(app.config['UPLOAD_PATH'], filename))
        if not validate(filepath):
            errmsg = "Invalid file format. Please choose a FASTA format file."
            os.remove(filepath)
            return render_template('index.html', errmsg=errmsg)
    return redirect(url_for('data',filename=filename))

@app.route('/uploads/<filename>')
def upload(filename):
    return send_from_directory(app.config['UPLOAD_PATH'], filename)

@app.route('/data/<filename>')
def data(filename):
    xmlfilename=parsefasta(filename)
    print(xmlfilename)
    #xmlfilename="./scripts/scriptxml_0.xml"  #For testing
    items=parsexml(xmlfilename)
    return render_template('data.html',items=items, filename=filename)

if __name__ == '__main__':
   app.run(debug = True)