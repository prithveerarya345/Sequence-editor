from flask import Flask, render_template, request, jsonify
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def index():
    result = None
    if request.method == "POST":
        seq = request.form.get("sequence")
        action = request.form.get("action")

        if action == "Complement Sequence":
            result = complement_sequence(seq)
        elif action == "Reverse Sequence":
            result = reverse_sequence(seq)
        elif action == "Translate Sequence":
            result = translate_sequence(seq)
        elif action == "Remove Numbers":
            result = rem_int(seq)
        elif action == "Remove Spaces":
            result = rem_sp(seq)
        elif action == "Remove linebreaks":
            result = rem_br(seq)
        else:
            result = seq  # Return the original sequence if no action is provided

        return jsonify(result=result)

    return render_template("index.html")

@app.route("/blast", methods=['GET', 'POST'])
def blast_page():
    if request.method == 'POST':
        sequence = request.form['sequence']
        blast_output = run_blast_search(sequence)
        return render_template('output.html', blast_output=blast_output)
    return render_template('blast.html')


def complement_sequence(seq):
    complement = ''
    for n in seq:
        if n.upper() == "A":
            complement += "T"

        if n.upper() == "T":
            complement += "A"

        if n.upper() == "G":
            complement += "C"

        if n.upper() == "C":
            complement += "G"

    return complement

def reverse_sequence(seq):
    rev = seq[::-1]
    return rev

def translate_sequence(seq):
    seq_obj = Seq(seq)
    return str(seq_obj.translate())

def rem_int(seq):
    res = ''
    for s in seq:
        if not s.isnumeric():
            res += s
        else:
            res = res

    return res

def rem_sp(seq):
    seq = seq.replace(" ", "")
    return seq

def rem_br(seq):
    seq = seq.replace("\n", "").replace("\r", "")
    return seq

def calculate_gc_content(seq):
    seq_obj = Seq(seq)
    gc_content = seq_obj.gc_content()
    return f"GC content: {gc_content:.2%}"

def run_blast_search(sequence):
    
    try:
        # Perform BLASTN search
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence, hitlist_size=1)

        # Parse the BLAST results
        blast_record = NCBIXML.read(result_handle)

        # Print the best matching hit
        if blast_record.alignments:
            best_alignment = blast_record.alignments[0]
            best_hsp = best_alignment.hsps[0]
            blast_output = (
                "=====Best Alignment=====\n"
                f"Sequence: {best_alignment.title}\n"
                f"Length: {best_alignment.length}\n"
                f"E-value: {best_hsp.expect}\n"
                f"Score: {best_hsp.score}\n"
                f"Query: {best_hsp.query}\n"
                f"Match: {best_hsp.match}\n"
                f"Sbjct: {best_hsp.sbjct}\n"
            )
            return blast_output
        else:
            return "No matching sequences found."

    except Exception as e:
        return "An error occurred:", str(e)
    
if __name__ == "__main__":
    app.run(host='127.0.0.1', port=8080, debug=True)
