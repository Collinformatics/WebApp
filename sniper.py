from flask import Flask, jsonify, render_template_string, request
import numpy as np
import pandas as pd
from figures import plotEntropy
import sys


app = Flask(__name__)

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
              'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def subsDefault():
    subs = {'VVLQAGTK': 17076,
            'LILQSVGA': 16783,
            'VALQSACW': 16331,
            'LNLQGILD': 15201,
            'IYLQALMP': 15010,
            'TSLQARKS': 14449,
            'AFLQAHFT': 13175,
            'VTLQCTYS': 13330,
            'VLLQAKQL': 12201,
            'LVLQANPC': 12010,
            'IMLQGVIW': 11449,
            'AGLQASAH': 11175,
            'MHLQSENE': 10330,
            'MVLQGDVN': 9201,
            'YGLQCNEV': 9010,
            'KCMQAQVQ': 8449,
            'GELQSWHF': 8175,
            'CDMQCMWG': 7330,
            'VWMQCSII': 7013,
            'VVFQCCNM': 2576,
            }

    return subs


def processData(subs, enzymeName, defaultSubs):
    subLen = len(next(iter(subs)))
    if defaultSubs:
        # Count: Substrates
        subsCounts = {}
        for sub in subs:
            if len(sub) == subLen:
                keepSub = True
                for aa in sub:
                    if aa not in AA:
                        keepSub = False
                        break
                if keepSub:
                    if sub in subsCounts.keys():
                        subsCounts[sub] += 1
                    else:
                        subsCounts[sub] = 1
    else:
        subsCounts = subs

    # Define: Substrate positions
    pos = [f'R{index + 1}' for index in range(subLen)]

    # Count AAs
    totalSubs = 0
    countedAA = pd.DataFrame(0, index=AA, columns=pos)
    for sub, count in subs.items():
        totalSubs += count
        for index, aa in enumerate(sub):
            countedAA.loc[aa, pos[index]] += count

    # Evaluate: AA probability
    probAA = countedAA / totalSubs

    # Figure: Entropy
    entropy, figEntropy = plotEntropy(probAA, AA, enzymeName)

    # Create dataset
    dataset = {}
    dataset['substrates'] = subsCounts
    dataset['probability'] = probAA
    dataset['entropy'] = entropy
    dataset['figEntropy'] = figEntropy

    return dataset



@app.route('/run', methods=['POST'])
def run():
    # Get file from the form
    try:
        substrates = request.files.get('textFile')
        loadFile = False
        if substrates:
            loadFile = True
            print('Subs: Load')

            # Read the contents of the file
            substrate = substrates.read().decode('utf-8')

            # Split the content by lines to get each substrate sequence
            substrates = substrate.splitlines()
        else:
            substrates = subsDefault()
            print('Subs: Default')

    except Exception as e:
        return jsonify({"error": f"Error A: {str(e)}"}), 400


    # Get other data from the form
    enzymeName = request.form.get('enzymeName')
    threshold = request.form.get('minS')
    NSubs = request.form.get('N')

    # Evaluate: Data
    dataset = processData(substrates, enzymeName, loadFile)

    result = {
        "fileReceived": substrates.filename if loadFile else False,
        "enzyme": enzymeName,
        "minS": threshold,
        "NSubs": NSubs,
        "figEntropy": dataset['figEntropy']
    }

    return jsonify(result)


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
                .container-description {
                    background-color: {{ grey }};
                    color: {{ white }};
                    font-size: 20px;
                    text-align: center;
                    border-radius: {{ borderRad }}px;
                    line-height: 1.4;
                    margin: 20px auto;
                    padding: 20px;
                }
                .container {
                    background-color: {{ grey }};
                    font-size: 20px;
                    text-align: center;
                    border-radius: {{ borderRad }}px;
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    padding: {{ padTB }}px {{ padSide }}px 
                             {{ marginB }}px {{ padSide }}px;
                    margin: {{ spacer }}px auto;
                    margin-bottom: {{ marginB }}px;
                }
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
                    border: 2px solid {{ greyDark }};
                    border-radius: {{ borderRad }}px;
                    width: 100%;
                    max-width: 350px;
                    padding: {{ padInput }}px;
                    margin-bottom: {{ padTB }}px;
                }
                button {
                    background-color: {{ green }};
                    color: {{ white }};
                    font-size: {{ fontSize }}px;
                    cursor: pointer;
                    border: none;
                    border-radius: {{ borderRad }}px;
                    padding: 8px 20px;
                    margin-top: {{ marginButton }}px;
                    margin-bottom: {{ marginButton }}px;
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
                }
                .form-group {
                    width: 100%;
                    max-width: 350px;
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
                    padding: {{ padInput }}px;
                    margin-bottom: {{ spacer }}px;
                }
            </style>
        </head>
        <script>
            function evaluateForm(event) {
                event?.preventDefault(); // Prevent any default form action
        
                const enzymeName = document.getElementById('enzymeName').value;
                const minS = document.getElementById('minS').value;
                const N = document.getElementById('N').value;
                const textFileInput = document.querySelector('input[name="textFile"]');
                const textFile = textFileInput.files[0];
        
                if (!enzymeName || !minS || !N) {
                    alert("Please fill out all required fields.");
                    return;
                }
        
                const formData = new FormData();
                formData.append('enzymeName', enzymeName);
                formData.append('minS', minS);
                formData.append('N', N);
                if (textFile) {
                    formData.append('textFile', textFile);
                }
        
                fetch('/run', {
                    method: 'POST',
                    body: formData
                })
                .then(response => {
                    if (!response.ok) throw new Error("Network response was not ok");
                    return response.json();
                })
                .then(data => {
                    const newDiv = document.createElement('div');
                    newDiv.className = 'container';
        
                    let headerContent = '';
                    if (data.fileReceived) {
                        headerContent = `
                            <p><strong>Enzyme:</strong> ${enzymeName}</p>
                            <p><strong>Min:</strong> ${minS}</p>
                            <p><strong>Top N:</strong> ${N}</p>
                        `;
                    } else {
                        headerContent = `<p>Use default substrates.</p>`;
                    }
        
                    newDiv.innerHTML = `
                        <div class="div-header">Results:</div>
                        ${headerContent}
                    `;
                    
                    // Check if figEntropy is available and add it as an image
                    if (data.figEntropy) {
                        const imgElement = document.createElement('img');
                        imgElement.src = 'data:image/png;base64,' + data.figEntropy;
                        imgElement.alt = 'Entropy Plot';
                        imgElement.style.maxWidth = '100%';  // Optional: control image size
                        newDiv.appendChild(imgElement);
                    }
            
                    document.body.appendChild(newDiv);
                })
                .catch(error => {
                    console.error('Error:', error);
                    const errorDiv = document.createElement('div');
                    errorDiv.className = 'container';
                    errorDiv.innerHTML = `<p style="color: red;">
                        Error fetching results: ${error.message}</p>`;
                    document.body.appendChild(errorDiv);
                });
            }
        </script>
        <body>
            <h1>{{ title|safe }}</h1>
            <h2>{{ header|safe }}</h2>
            <div class="container-description">
                <p>{{ ct1|safe }}</p>
                <p>{{ ct2|safe }}</p>
                <p>{{ equation|safe }}</p>
                <p>{{ ct3|safe }}</p>
                <p>{{ ct4|safe }}</p>
            </div>
            <div class="container">
                <form action="/run" method="POST" 
                      enctype="multipart/form-data" class="input-form">
                    <div class="div-header">Upload Text File:</div>
                    <input type="file" name="textFile" accept=".txt"
                    class="input-form"> <!-- removed 'required' -->
                    
                    <div class="div-header">Experiment Parameters:</div>
                    
                    <div class="form-group">
                        <label for="enzymeName">Enzyme Name:</label>
                        <input type="text" id="enzymeName" name="enzymeName" required
                            value=name>
                    </div>
                    
                    <div class="form-group">
                        <label for="minS">Entropy Threshold:</label>
                        <input type="number" id="minS" name="minS" 
                        step="0.1" value="0.6" required>
                    </div>
                    
                    <div class="form-group">
                        <label for="N">Number of Substrates:</label>
                        <input type="number" id="N" name="N" value="30" required>
                    </div>
                    
                    <button type="button" onclick="evaluateForm()">Submit</button>
                </form>
            </div>
        </body>
    </html>
    ''',
    black='#151515', white='FFFFFF',
    grey='#303030', greyDark='#202020',
    green='#23FF55', greenLight='#1DE64A',

    spacer=20, spacerMini=5,
    padSide=50, padTB=30, padInput=8,
    marginB=12, marginButton=12,
    fontSize=18,
    borderRad=5,


    title="Specificity Network Identification via Positional Entropy based Refinement "
          "(SNIPER)",

    header="Modeling Enzyme Specificity",

    ct1="This program will take substrates for a given enzyme and identify the "
        "Motif, of the recognition sequence within the larger protein sequence. "
        "The Motif is identified by the positions in the substrate that have "
        "Entropy scores (∆S) that exceed a minS value.",

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
