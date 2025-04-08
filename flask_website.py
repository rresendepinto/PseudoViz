# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 15:39:29 2022

@author: Rui
"""

from urllib import request
from flask import Flask, flash, redirect, url_for, render_template, request, session
import os
from random import randint
import json
from Bio import SeqIO, Seq
from werkzeug.utils import secure_filename
from utils import *
import dill
from math import floor
from datetime import timedelta
import secrets
from uuid import uuid4
import time
import subprocess as sp
import shutil
from flask_jsonpify import jsonify

app = Flask(__name__)

app.secret_key = secrets.token_bytes(32)
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days = 7)

app.static_folder = 'static'

UPLOAD_FOLD = 'UPLOAD_FOLDER'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLD

cds_dict = {}

@app.route("/")
def home():
    
        return render_template("index.html")

@app.route("/testd3")
def test_d3():
    
        return render_template("test_d3.html")

@app.route("/uploader", methods=["GET", "POST"])
def uploader():
    global analysis_name
    global session_folder
  
  

    if request.method == 'POST':
        f = request.files['file']
        
        analysis_name = request.form['analysis_name']

        session['number']= str(uuid4())
        session_folder = f'{analysis_name}_{session["number"]}'
        session['session_folder'] = session_folder
        os.mkdir(UPLOAD_FOLD + '/' + session_folder)
        with open(UPLOAD_FOLD +  '/' + session_folder + '/test_upload_fold', 'wt') as handle:
            handle.write('Hello World!')
        #os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], f'{session_folder}'))

        f.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename( 'uploaded_reference_file.fa')))
        exons_data = request.files['exons_data']
        exons_data.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename( 'exons_data.json')))
        genomic_file = request.files['genomic_fasta']
        genomic_file.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename('genomic_fasta.fa')))
  
        return redirect(url_for("stats"))
    else:
        return "No file uploaded"
        
@app.route("/uploader_macse", methods=["GET", "POST"])
def uploader_macse():
       
    macse_analysis_list[session['session_folder']] = []

    if request.method == 'POST':
        macse_analysis_name = request.form['macse_analysis_name']
        fasta_aa = request.files['macseanalysis_aa_ordered']
        fasta_aa.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename(f'{macse_analysis_name}_macseanalysis_aa_ordered_file.fa')))
        fasta_nt = request.files['macseanalysis_nt_ordered']
        fasta_nt.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename(f'{macse_analysis_name}_macseanalysis_nt_ordered_file.fa')))
        #config_file = request.files['config']
        #config_file.save(secure_filename("config_params.json"))
        exon_data_dill = request.files['first_step_data']
        exon_data_dill.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename(f'first_step_data_dill.pkl')))
        macse_analysis_list[session['session_folder']].append(macse_analysis_name)
  
        return redirect(url_for("macse_alignments"))
    else:
        return "No file uploaded"

@app.route("/upload_example", methods=["GET", "POST"])
def upload_example():
    global analysis_name
    global session_folder
    

    if request.method == 'POST':
        
        
        analysis_name = 'Example Data'

        session['number']= str(uuid4())
        session_folder = f'{analysis_name}_{session["number"]}'
        session['session_folder'] = session_folder
        os.mkdir(UPLOAD_FOLD + '/' + session_folder)

        macse_analysis_list[session['session_folder']] = []
        macse_analysis_list[session['session_folder']].append("Example")

        directory = os.path.dirname(os.path.abspath(__file__))


        shutil.copy(os.path.join(directory, 'example_data/CYP2J19_reference.fa'), 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'uploaded_reference_file.fa')
        shutil.copy(os.path.join(directory, 'example_data/exon_alns.json'), 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/'+ 'exons_data.json')
        shutil.copy(os.path.join(directory, 'example_data/CYP2J19_genomic_regions.fna'), 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'genomic_fasta.fa')
        shutil.copy(os.path.join(directory, 'example_data/MACSE_test/macseanalysis_aa_ordered.fasta'), f'UPLOAD_FOLDER/{session["session_folder"]}/Example_macseanalysis_aa_ordered_file.fa')
        shutil.copy(os.path.join(directory, 'example_data/MACSE_test/macseanalysis_nt_ordered.fasta'), f'UPLOAD_FOLDER/{session["session_folder"]}/Example_macseanalysis_nt_ordered_file.fa')
        shutil.copy(os.path.join(directory, 'example_data/MACSE_test/first_step_data_dill.pkl'), 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'first_step_data_dill.pkl')
  
        return redirect(url_for("stats"))
    else:
        return "No file uploaded"

@app.route("/submit_macse_output")
def submit_macse_output():
    pass
    #return render_template()
    

@app.route("/exon_alignment_display", methods=['GET', 'POST'])
def exon_alignment_display():
    exon_n = '1'

    try:

        with open( 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'exons_data.json', "r") as handle:
            y = json.load(handle)

        handle = open( 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'uploaded_reference_file.fa', 'r')
        list_exs_cds = list(SeqIO.parse(handle, 'fasta'))

    except:

        return render_template('MissingFile.html')

    global numtotalexons
    numtotalexons = len(list_exs_cds) - 1

    
    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)
    
    pos_exons = get_pos_exons(list_exs_cds, cds)

    len_exon = pos_exons['Exon_' + exon_n][1] - pos_exons['Exon_' + exon_n][0] + 1
    
    my_json = json.dumps(global_pstops)

    #amino_dict = get_amino_from_cds_dict(cds_dict)
    

    def get_global_pos(exon_n, exon_pos) :
        #convert exon pos to global pos
        start_exon = pos_exons['Exon_' + str(exon_n)][0]

        return exon_pos + start_exon

    #throw error message in case analysis name is none

    if request.method == 'POST': 
        if request.form.get('options'):
            
            exon_n = request.form.get('options')
            len_exon = pos_exons['Exon_' + exon_n][1] - pos_exons['Exon_' + exon_n][0] + 1

        
            return render_template("exon_alignment_display.html", targets=y, len=len_exon, upper_str=upper_str,
                                     numtotalexons=numtotalexons, exon_n=exon_n, pos_exon=pos_exons['Exon_' + exon_n][0], analysis_name=analysis_name,
                                      global_frameshifts=global_frameshifts, global_pstops=global_pstops, get_global_pos=get_global_pos)



    elif request.method == 'GET':
        
        return render_template("exon_alignment_display.html", targets=y, len=len_exon, upper_str=upper_str, 
                                numtotalexons=numtotalexons, exon_n=exon_n, pos_exon=pos_exons['Exon_' + exon_n][0], 
                                analysis_name=analysis_name, global_frameshifts=global_frameshifts, 
                                global_pstops=global_pstops, get_global_pos=get_global_pos)

@app.route("/exon_painter/<name>")
def exon_painter(name=None):
    
    with open( 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'genomic_fasta.fa', "r") as handle:
        genomic_fasta = {}
        #this will deal with duplicate sequences in fasta
        records= list(SeqIO.parse(handle, 'fasta'))
        for seq in records:
            if seq.id not in genomic_fasta.keys():
                genomic_fasta[seq.id] = seq


    genomic_region = genomic_fasta[name].seq

    with open( 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'exons_data.json', "r") as handle:
        y = json.load(handle)

    target = y[name]


    target_exons_pos = []
    correct_splice_sites = []
    wrong_splice_sites = []

    for x in target.keys():

        target_exons_pos.append((target[x][7], target[x][8]))
        
        if x != "1":
            if target[x][3].upper() == 'AG':

                correct_splice_sites.append(target[x][7]-2)
                correct_splice_sites.append(target[x][7]-1)

            else: 

                wrong_splice_sites.append(target[x][7]-2)
                wrong_splice_sites.append(target[x][7]-1)

        if x != str(numtotalexons):
            if target[x][4].upper() == 'GT' or target[x][4].upper() == 'GC':
                
                correct_splice_sites.append(target[x][8]+1) 
                correct_splice_sites.append(target[x][8]+2)

            else:

                wrong_splice_sites.append(target[x][8]+1) 
                wrong_splice_sites.append(target[x][8]+2)



    def check_if_pos_is_exon(pos):

        pos_is_exon = False

        for x in target_exons_pos:

            if pos >= x[0] and pos <= x[1]:

                pos_is_exon = True

        return pos_is_exon

    def is_correct_splice_site(pos):
        
        if pos in correct_splice_sites:

            return True
        else:
            return False
        
    def is_wrong_splice_site(pos):
        
        if pos in wrong_splice_sites:

            return True
        else:
            return False

    len_region=len(genomic_region)

    number_of_rows= len_region//150
    last_row = len_region % 150

    list_start_end = []

    cur_start = 0
    for i in range(number_of_rows + 1):
        list_start_end.append((cur_start, cur_start + 150))
        cur_start = cur_start + 150


    return render_template('exon_painter.html', genomic_region=genomic_region,  number_of_rows= number_of_rows, last_row = last_row, list_start_end=list_start_end,
                                    is_pos_exon = check_if_pos_is_exon, is_correct_splice_site = is_correct_splice_site, is_wrong_splice_site=is_wrong_splice_site)   

@app.route("/stats")
def stats():

    
    table_frameshifts = pd.DataFrame()
    table_premature_stops = pd.DataFrame()

    try:

        handle = open('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'uploaded_reference_file.fa', 'r')
        list_exs_cds = list(SeqIO.parse(handle, 'fasta'))
        

        with open('UPLOAD_FOLDER' + '/' + session['session_folder'] +  '/' + 'exons_data.json', "r") as handle:
            exon_data = json.load(handle)
    
    except:

        return render_template('MissingFile.html')
    
    pseudoindex_list = []



    for sp_id in exon_data.keys():
        try:
            sp_id_table_frameshifts = create_mutations_table(exon_data[sp_id], sp_id)

            sp_id_table_frameshifts2 = sp_id_table_frameshifts.copy()

            sp_id_table_pstops, exonscat_withgaps = create_pstop_table(exon_data[sp_id], list_exs_cds, sp_id, sp_id_table_frameshifts)
            
            cds_dict[sp_id] = exonscat_withgaps

            

            if table_frameshifts.empty:
                table_frameshifts = sp_id_table_frameshifts
                
            else:
                table_frameshifts  = pd.concat([table_frameshifts, sp_id_table_frameshifts])

            if table_premature_stops.empty:
                table_premature_stops = sp_id_table_pstops

            else:
                table_premature_stops  = pd.concat([table_premature_stops, sp_id_table_pstops])

            
            list_for_table = pseudoindex_stats(sp_id_table_frameshifts, sp_id_table_pstops, exon_alns=exon_data[sp_id], list_exs_cds = list_exs_cds)
            list_for_table.insert(0, sp_id)
            pseudoindex_list.append(list_for_table)
            
            #I need to iterate through the pstops and frameshifts tables to alter the positions of the premature stops
                #cause if there is an insertion before the pstop in the same exon, the wrong position will be highlighed as a pstop

            
            
            #create list from dataframe
            global_pstops[sp_id] = create_pstops_list(sp_id_table_pstops, list_exs_cds, cds= str(list_exs_cds[len(list_exs_cds) - 1].seq))
            #print(global_pstops[sp_id])
            
            new_global_pstops_sp_id = []
            for global_pstop in global_pstops[sp_id]:
                new_global_pstop = global_pstop
                exon_pstop = int(global_pstop.split('Exon_')[1].split('_')[0])
                pos_pstop_on_exon =  int(global_pstop.split('Exon_')[1].split('_')[1])
                #print(exon_pstop, pos_pstop_on_exon)
                
                for index, row in sp_id_table_frameshifts2.iterrows():
                    if row['type'] == 'insertion' and int(row['exon']) == exon_pstop and int(row['start']) < pos_pstop_on_exon:
                        print('removing:' + global_pstop)
                        new_global_pstop = '_'.join(global_pstop.split('_')[:-1]) + '_' + str(pos_pstop_on_exon + row['len'])
                    
                
                new_global_pstops_sp_id.append(new_global_pstop)
            global_pstops[sp_id] = new_global_pstops_sp_id
            
        except:
            print(sp_id)
        
    
    
    pseudoindex_table = pd.DataFrame(pseudoindex_list, columns=['target', 'pseudoindex', 'shifted_frame%', 'truncated_frame%', 'missing_frame%']).round(2)
    
    general_stats_dict = general_stats(exon_data, list_exs_cds)
    general_stats_table = pd.DataFrame.from_dict(general_stats_dict).T.reset_index()
    

    general_stats_columns = ['Target', "Avg. Exon Align. Pid", 'Avg. Exon Size', 'Min Exon Size', 'Larg. Exon size', 'Splice site integrity', 'Aligning exons']


    return render_template("Stats.html", tables=[general_stats_table.to_html(classes='data'), 
                            table_frameshifts.to_html(classes='data'), pseudoindex_table.to_html(classes='data')], 
                            titles_frameshifts=table_frameshifts.columns.values, row_data_frameshifts=list(table_frameshifts.values.tolist()), 
                            titles_pstops = table_premature_stops.columns.values, row_data_pstops = list(table_premature_stops.values.tolist()),
                            titles_general_stats = general_stats_columns, row_data_general_stats= list(general_stats_table.values.tolist()), 
                            titles_pseudoindex = pseudoindex_table.columns.values, row_data_pseudoindex = list(pseudoindex_table.values.tolist()),
                            zip=zip)
        



    #return render_template("stats.html")

@app.route("/processing", methods=["GET", "POST"])
def processing():

    global analysis_name
    global session_folder
  
  

    if request.method == 'POST':
        f = request.files['file']
        
        analysis_name = request.form['analysis_name']

        session['number']= str(uuid4())
        session_folder = f'{analysis_name}_{session["number"]}'
        session['session_folder'] = session_folder
        os.mkdir(UPLOAD_FOLD + '/' + session_folder)

        f.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename( 'uploaded_reference_file.fa')))
        genomic_file = request.files['genomic_fasta']
        genomic_file.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename('genomic_fasta.fa')))
        main_path_name = 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/'
        reference = main_path_name + 'uploaded_reference_file.fa'
        genome_regions = main_path_name  + 'genomic_fasta.fa'
        pseudochecker_process = f'python pseudochecker2.0/pseudochecker.py --file_exs_cds  {reference} --file_genomic {genome_regions}  --analysis_name {str(analysis_name)} --main_path {main_path_name} --skip_MACSE'
        process1 = subprocess.Popen(pseudochecker_process, stdout=subprocess.PIPE, shell=True)
        output, error = process1.communicate()
        process1.wait()
        print(pseudochecker_process)
        print(output)
        print(error)
        while os.path.exists(main_path_name + 'exon_alns.json') == False:
            if os.path.exists(main_path_name + 'stopped.txt'):
                #add here an error information page
                break
            time.sleep(5)

        shutil.copyfile(main_path_name + 'exon_alns.json', 'UPLOAD_FOLDER' + '/' + session['session_folder'] + '/' + 'exons_data.json')
        #exons_data.save(os.path.join('UPLOAD_FOLDER' + '/' + session['session_folder'] + '/', secure_filename( 'exons_data.json')))
        
  
        return redirect(url_for("stats"))
    else:
        return "No file uploaded"

   

    return redirect(url_for("upload_example"))

@app.route("/macse_alignments")
def macse_alignments():
    global analysis_name

    if session['session_folder'] not in macse_analysis_list.keys():

         return render_template("submit_MACSE_alignments.html")
    
    else:

        if len(macse_analysis_list[session['session_folder']]) == 0:

            return render_template("submit_MACSE_alignments.html")

    
        macse_analysis_name = macse_analysis_list[session['session_folder']][-1]
        
        results_aa = open( f'UPLOAD_FOLDER/{session["session_folder"]}/{macse_analysis_name}_macseanalysis_aa_ordered_file.fa', 'r')
        results_aa = list(SeqIO.parse(results_aa, 'fasta'))  # Parses the MACSE Alignment in AminoAcids


        

        results_nt = open(f'UPLOAD_FOLDER/{session["session_folder"]}/{macse_analysis_name}_macseanalysis_nt_ordered_file.fa', 'r')
        results_nt = list(SeqIO.parse(results_nt, 'fasta'))  # Parses the MACSE Alignment in Nucleotides
        

        

        with open(f'UPLOAD_FOLDER/{session["session_folder"]}/first_step_data_dill.pkl', 'rb') as f:
            data_needed_for_macse =  dill.load(f)

       

        cdsposglobal = data_needed_for_macse[0]
        partialsequences = data_needed_for_macse[1]
        excludedsequences = data_needed_for_macse[2]
        no_additional_cds = data_needed_for_macse[3]
        numtotalexons = data_needed_for_macse[4]
        exonoccupancy  = data_needed_for_macse[5]
        additionalcds = data_needed_for_macse[6]
        nrglobal = data_needed_for_macse[7]

        excludedsequences = str(excludedsequences).replace(']', '')
        excludedsequences = excludedsequences.replace('[', '')

        partialsequences_str = []
        for i in partialsequences:
            if i not in excludedsequences:
                partialsequences_str.append(i)
        partialsequences_str = str(partialsequences_str).replace(']', '')
        partialsequences_str = partialsequences_str.replace('[', '')

        aln_len = alignmentlength(results_aa)
        noseqs = nosequences(results_aa)
        avg_pairwise_ident = averagepairwiseaaidentity(results_aa)
        aa_ident_sites = aaidenticalsites(results_aa)

       
        
        def get_color(j):
            colorlist=['#CCFFC4', '#DBEBFF', '#E3E3E3']
            color_idx = floor(((j % 9) - 1) / 3)
        
            return colorlist[color_idx]

        return render_template("MACSE_alignments.html", range=range, len=len, str=str, results_nt=results_nt, results_aa = results_aa, excludedsequences= excludedsequences, 
                    no_additional_cds=str(no_additional_cds), aln_len=aln_len, noseqs = noseqs, avg_pairwise_ident = avg_pairwise_ident, aa_ident_sites=aa_ident_sites, analysis_name=macse_analysis_name,
                    replace=replace, int=int, get_color=get_color)

ALLOWED_EXTENSIONS = {'json'}

@app.route("/dendrogram", methods=['GET'])
def dendrogram():
    print("Dendrogram page route was called")
    return render_template("dendrogram.html")

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/upload_tree', methods=['POST'])
def upload_tree():
    if 'datafile' not in request.files:
        return {'error': 'No file part'}, 400
    file = request.files['datafile']
    
    if file.filename == '':
        return {'error': 'No selected file'}, 400
    if file and allowed_file(file.filename):
        
        data = json.load(file)  # Load JSON data from the file
        
        return jsonify(data)    # Send JSON data back to the client
    else:
        return {'error': 'File type not allowed'}, 400

if __name__ == "__main__":

    seq1=""
    seq2=""

    
    
    list_exs_cds = []
    y=None
    
    numtotalexons = len(list_exs_cds) - 1

   
    analysis_name = None
    macse_analysis_list  = {}
    session_folder = None

    global_frameshifts = {}
    global_pstops = {}

    app.run(debug=True, use_reloader=True, host="0.0.0.0")
    #app.run(debug=False, host="0.0.0.0")
    #app.run(debug=False,  port='8088', host='193.136.51.225')

