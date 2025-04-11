from flask import Flask, render_template_string


inTitle = ('Specificity Network Identification via '
           'Positional Entropy based Refinement (SNIPER)')

app = Flask(__name__)


@app.route('/upload', methods=['POST'])
def upload():
    fastaFile = request.files['fastaFile']
    # process the uploaded file...
    return f"Uploaded: {fastaFile.filename}"


@app.route('/')
def home():
    return render_template_string('''
    <html>
        <head>
            <title>SNIPER</title>
            <style>
                body {
                    background-color: {{black}};
                    color: #f0f0f0;
                    font-family: Arial, sans-serif;
                    padding: 40px;
                }
                h1 {
                    text-align: center;
                    font-size: 35px;
                    padding: 0px;
                }
                h2 {
                    text-align: center;
                    font-size: 25px;
                    padding: 10px;
                }
                p {
                    font-size: 18px;
                    text-align: center;
                    padding-top: 5px;
                    padding-bottom: 0px;
                }
                .description-container {
                    background-color: {{grey}};
                    color: {{white}};
                    font-size: 20px;
                    margin-top: 30px;
                    border-radius: 5px;
                    text-align: center;
                    padding: 5px 20px;
                    text-align: center;
                    margin: 20px auto;
                    line-height: 1.2;
                }
                .upload-container {
                    background-color: {{grey}};
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    border-radius: 10px;
                    width: 420px;
                    margin: 0px auto;
                    box-shadow: 0 0 10px rgba(0, 0, 0, 0.5);
                    padding: 30px 50px 30px 50px;
                }
                .div-container {
                    background-color: {{grey}};
                    color: {{white}};
                    font-size: 20px;
                    width: 420px; /* Added 'px' */
                    margin: 30px auto; /* Centers it horizontally */
                    border-radius: 5px;
                    text-align: center;
                    padding: 20px;
                }        
                .upload-header {
                    color: #23FF55;
                    font-size: 20px;
                    font-weight: normal;
                    padding-top: 0px;
                    padding-bottom: 15px;
                }
                input[type="file"] {
                    background-color: {{black}};
                    color: #101010;
                    font-size: 18px;
                    margin-bottom: 20px;
                    border-radius: 6px;
                    width: 100%;
                    max-width: 350px;
                    border: 2px solid {{greyDark}};
                    padding: 8px;
                }
                .input-form {
                    display: flex;
                    flex-direction: column;
                    align-items: center;  /* centers the inputs */
                }
                .input-form input {
                    background-color: {{ black }};
                    color: #f0f0f0;
                    font-size: 18px;
                    text-align: center;
                    margin-bottom: 20px;
                    border-radius: 6px;
                    width: 100%;
                    max-width: 350px;
                    border: 2px solid #303030;
                    padding: 10px;
                }
                button {
                    background-color: {{green}};
                    color: {{white}};
                    font-size: 18px;
                    cursor: pointer;
                    border: none;
                    border-radius: 5px;
                    padding: 8px 20px;
                }
                button:hover {
                    background-color: {{greenLight}};
                }
            </style>
        </head>
        <body>
            <h1>{{ title|safe }}</h1>
            <h2>{{ header|safe }}</h2>
            <div class="description-container">
                <p>{{ ct1|safe }}</p>
                <p>{{ ct2|safe }}</p>
                <p>{{ equation|safe }}</p>
                <p>{{ ct3|safe }}</p>
                <p>{{ ct4|safe }}</p>
            </div>
            <form class="upload-container" method="POST" 
                    action="/upload" enctype="multipart/form-data">
                <div class="upload-header">Upload FASTA File:</div>
                <input type="file" name="fastaFile" accept=".fasta,.fa" required>
                <button type="submit">Upload</button>
            </form>
            <div class="div-container">
            
            
                <form action="/run" method="POST" class="input-form">
                    <label for="enzymeName">Enzyme Name:</label><br>
                    <input type="text" id="enzymeName" name="enzymeName" required><br><br>

                    <label for="threshold">Entropy Threshold:</label><br>
                    <input type="number" id="threshold" name="threshold" step="0.1" required><br><br>
            
                    <label for="topN">Top N Substrates:</label><br>
                    <input type="number" id="topN" name="topN" required><br><br>
            
                    <button type="submit">Evaluate</button>
                </form>
                
                
            </div>
        </body>
    </html>
    ''',
    black='#151515', white='FFFFFF',
    grey='#303030', greyDark='#202020',
    green='#23FF55', greenLight='#25EC58',

    title=inTitle,

    header="Modeling Enzyme Specificity",

    ct1="This program will take substrates for a given enzyme and identify the "
        "Motif, of the recognition sequence within the larger protein sequence. "
        "The Motif is identified by the positions in the substrate that have "
        "Entropy scores (∆S) that exceed a threshold value.",

    ct2="∆S is evaluated at each position in the substrate sequence and is found by "
        "the difference between the Maximum Entropy (S<sub>Max</sub>) and the "
        "Shannon Entropy (S<sub>Shannon</sub>)",

    equation="∆S = S<sub>Max</sub> - S<sub>Shannon</sub> = log<sub>2</sub>(20) - "
            "∑(-prob<sub>AA</sub> * log<sub>2</sub>(prob<sub>AA</sub>))",

    ct3="Once we have identified the motif, we can select the most common AA sequences "
        "and use them as input for a Suffix Tree. This will plot the preferred residues "
        "in descending order of Specificity, as determined by the ∆S values.",
    ct4="The tree will reveal the Specificity Network of your enzyme. "
        "In other words, we will be able to identify the AA preferences at a given "
        "location in the motif when a specific AA is at an other position."
)



if __name__ == '__main__':
    app.run(debug=True)
