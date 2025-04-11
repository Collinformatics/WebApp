from flask import Flask, render_template_string


inTitle = ('Specificity Network Identification via '
           'Positional Entropy based Refinement (SNIPER)')

# width: 420px;

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
                    background-color: {{ black }};
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
                    font-size: {{ fontSize }}px;
                    text-align: center;
                    padding-top: 5px;
                    padding-bottom: 0px;
                }
                .description-container {
                    background-color: {{ grey }};
                    color: {{ white }};
                    font-size: 20px;
                    text-align: center;
                    border-radius: {{ borderRad }}px;
                    line-height: 1.2;
                    margin: 20px auto;
                    padding: 15px 20px;
                }


                .container {
                    background-color: {{ grey }};
                    font-size: 20px;
                    text-align: center;
                    border-radius: {{ borderRad }}px;
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    box-shadow: 0 0 10px rgba(0, 0, 0, 0.5);
                    margin: {{ spacer }}px auto;
                    margin-bottom: {{ marginB }}px;
                    padding: {{ 30 }}px 50px 30px 50px;
                }

                .div-container {
                    background-color: {{ grey }};
                    color: {{ white }};
                    font-size: 20px;
                    text-align: center;
                    border-radius: {{ borderRad }}px;
                    margin: {{ spacer }}px auto;
                    margin-bottom: {{ marginB }}px;
                    padding: {{ 30 }}px 50px 30px 50px;
                } /* width: 420px; */


                .div-header {
                    color: #23FF55;
                    font-size: 20px;
                    font-weight: normal;
                    padding-top: 0px;
                    padding-bottom: 15px;
                }
                input[type="file"] {
                    background-color: {{ black }};
                    color: #101010;
                    font-size: {{ fontSize }}px;
                    margin-bottom: {{ marginB }}px;
                    border: 2px solid {{ greyDark }};
                    border-radius: {{ borderRad }}px;
                    width: 100%;
                    max-width: 350px;
                    padding: 8px;
                }
                button {
                    background-color: {{ green }};
                    color: {{ white }};
                    font-size: {{ fontSize }}px;
                    cursor: pointer;
                    border: none;
                    border-radius: {{ borderRad }}px;
                    padding: 8px 20px;
                }
                button:hover {
                    background-color: {{ greenLight }};
                }
                .input-form {
                    align-items: center;
                    text-align: center; /* center form text */
                    display: flex;
                    flex-direction: column;
                    font-size: {{ fontSize }}px;
                    margin-bottom: 12px;
                }
                .form-group {
                    width: 100%;
                    max-width: 350px;
                    margin-bottom: 12px;
                }
                .form-group label {
                    display: block;
                    margin-bottom: {{ spacerMini }}px;
                    font-size: {{ fontSize }}px;
                    text-align: center;
                }
                .form-group input {
                    font-size: {{ fontSize }}px;
                    border-radius: {{ borderRad }}px;
                    border: 2px solid {{ greyDark }};
                    text-align: center;
                    width: 100%;
                    margin-bottom: {{ spacer }}px;
                    padding: 8px;
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
            <form class="container" method="POST" 
                    action="/upload" enctype="multipart/form-data">
                <div class="div-header">Upload FASTA File:</div>
                <input type="file" name="fastaFile" accept=".fasta,.fa" required>
                <button type="submit">Upload</button>
            </form>
            <div class="container">
                <form action="/run" method="POST" class="input-form">
                    <div class="div-header">Experiment Parameters:</div>
                    <div class="form-group">
                        <label for="enzymeName">Enzyme Name:</label>
                        <input type="text" id="enzymeName" name="enzymeName" required>
                    </div>

                    <div class="form-group">
                        <label for="threshold">Entropy Threshold:</label>
                        <input type="number" id="threshold" name="threshold" 
                        step="0.1" required>
                    </div>

                    <div class="form-group">
                        <label for="topN">Select Number of Substrates:</label>
                        <input type="number" id="topN" name="topN" required>
                    </div>
            
                    <button type="submit">Evaluate</button>
                </form>
</div>
        </body>
    </html>
    ''',
    black='#151515', white='FFFFFF',
    grey='#303030', greyDark='#202020',
    green='#23FF55', greenLight='#1DE64A',

    spacer=20, spacerMini=5,
    padSide=20, padTB=20,
    marginB=12,
    fontSize=18,
    borderRad=5,


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
