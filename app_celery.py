from ast import parse
from fileinput import close
import os
from flask import Flask, render_template, request,redirect, url_for, Response,send_from_directory,jsonify
from celery import Celery
from werkzeug.utils import secure_filename
from werkzeug.datastructures import  FileStorage
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time

app = Flask(__name__)   
app.config['UPLOAD_PATH'] = 'uploads'

app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
# Celery configuration
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)
# Initialize Celery
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

def validate(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # Return True if this is fasta file format

def parsexml(filenames):
    files=[]
    items=[]
    i=0
    for filename in filenames:
        i+=1
        for record in NCBIXML.parse(open(filename)):
            if record.alignments:
                for align in record.alignments:
                    for hsp in align.hsps: 
                        item=dict({"query_id":record.query_id,"query_length":record.query_length,"sequence":align.title,"evalue":hsp.expect })
                        items.append(item)
        files.append(items)
    return i, files

@celery.task(bind=True)
def parsefasta(filename):
    i=0
    xmlfiles=[]
    filepath = os.path.join("./",app.config['UPLOAD_PATH'], filename) 
    with open(filepath) as fh:
             for seq_record in SimpleFastaParser(fh):
                i+=1
                xmlfilename = "results_"+str(i)+".xml"
                with open(xmlfilename, "w") as save_file:
                    print(seq_record[1])
                    result_handle = NCBIWWW.qblast("blastn", "nt", sequence=seq_record[1], format_type="XML") 
                    save_file.write(result_handle.read())
                xmlfiles.append(xmlfilename)
    return xmlfiles

@app.route('/longtask', methods=['POST'])
def longtask():
    task = parsefasta.apply_async()
    return jsonify({}), 202, {'Location': url_for('taskstatus',
                                                  task_id=task.id)}

@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = parsefasta.AsyncResult(task_id)
    if task.state == 'PENDING':
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)


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
    xmlfiles=parsefasta(filename)
    i,files=parsexml(xmlfiles)
    return render_template('data.html',files=files,filename=filename,i=i)

@app.route('/progress')
def progress():
	def generate():
		x = 0
		
		while x <= 100:
			yield "data:" + str(x) + "\n\n"
			x = x + 10
			time.sleep(0.5)

	return Response(generate(), mimetype= 'text/event-stream')


if __name__ == '__main__':
   app.run(debug = True)