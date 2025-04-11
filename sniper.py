from flask import Flask, jsonify, render_template_string, request
import pandas as pd
import sys


app = Flask(__name__)

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
              'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def subsDefault():
    subs = {'VVLQAGCK': 22076,
            'LILQSVGH': 21783,
            'VALQSGCK': 20331,
            'LNLQGGCK': 19201,
            'IGLQAGCK': 16310,
            'VGLQAGCK': 15449,
            'LFLQAGCK': 12375,
            'VVLQCGCK': 10330,
            'FVLQAAGE': 10330,
            'VVMQCACK': 2576,
            }
    return subs


def countSubstrates(subs, defaultSubs):
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

    # Evaluate: AA distributions
    probAA = evalAA(subs, subLen)

    # Create dataset
    dataset = {}
    dataset['substrates'] = subsCounts
    dataset['probability'] = probAA

    return dataset


def evalAA(subs, subLen):
    # Define: Substrate positions
    pos = [f'R{index+1}' for index in range(subLen)]

    # Count AAs
    totalSubs = 0
    countedAA = pd.DataFrame(0, index=AA, columns=pos)
    for sub, count in subs.items():
        totalSubs += count
        for index, aa in enumerate(sub):
            countedAA.loc[aa, pos[index]] += count

    # Evaluate: AA probability
    probAA = countedAA / totalSubs

    # Evaluate: Entropy
    

    return probAA


@app.route('/run', methods=['POST'])
def run():
    # Get file from the form
    try:
        substrates = request.files['textFile']
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

    dataset = countSubstrates(substrates, loadFile)

    # Get other data from the form
    enzymeName = request.form.get('enzymeName')
    threshold = request.form.get('minS')
    NSubs = request.form.get('N')

    result = {
        "fileReceived": substrates.filename if loadFile else False,
        "substrates": substrates,
        "enzyme": enzymeName,
        "minS": threshold,
        "NSubs": NSubs,
    }

    print('Upload data')
    return jsonify(result)


@app.route('/')
def home():
    return render_template_string('''
    <html>
        <head>
            <title>SNIPER</title>
            <script>
                function evaluateForm() {
                    const enzymeName = document.getElementById('enzymeName').value;
                    const minS = document.getElementById('minS').value;
                    const N = document.getElementById('N').value;
                    const substrates = document.getElementById('substrates').value;
                    document.querySelector('input[name="textFile"]').files[0];
                
                    const formData = new FormData();
                    formData.append('enzymeName', enzymeName);
                    formData.append('minS', minS);
                    formData.append('N', N);
                    const textFile = document.querySelector('input[name="textFile"]').files[0];
                    formData.append('textFile', textFile);
                
                    fetch('/run', {
                        method: 'POST',
                        body: formData  // Don't set Content-Type manually
                                        // browser will handle it
                    })
                    .then(response => response.json())
                    .then(data => {
                        const newDiv = document.createElement('div');
                        newDiv.className = 'container';
                        const formattedJSON = JSON.stringify(data, null, 2);
                    
                        let headerContent = '';
                        if (data.fileReceived) {
                            headerContent = `
                                <p><strong>Enzyme:</strong> ${enzymeName}</p>
                                <p><strong>Min:</strong> ${minS}</p>
                                <p><strong>Top N:</strong> ${N}</p>
                            `;
                    
                            // If substrates are counted in the response
                            contentToShow = data.substrates || {};
                        } else {
                            headerContent = `
                                <p>No file received.</p>
                            substrates = subsDefault();
                            `;
                        }
                        const formattedJSON = JSON.stringify(contentToShow, null, 2);
                    
                        newDiv.innerHTML = `
                            <div class="div-header">Results:</div>
                            ${headerContent}
                            <pre style="text-align: left; 
                                 background-color: {{ black }}; 
                                 padding: 10px; 
                                 border-radius: {{ borderRad }}}px; 
                            </pre>
                        `;
                        document.body.appendChild(newDiv);
                    })

                    .catch(error => {
                        console.error('Error:', error);
                        const errorDiv = document.createElement('div');
                        errorDiv.className = 'container';
                        errorDiv.innerHTML = `<p style="color: red;">
                            Error fetching results: ${error.message || error}</p>`;
                        document.body.appendChild(errorDiv);
                    });
                }
            </script>
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
