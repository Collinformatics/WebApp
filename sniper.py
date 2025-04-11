from flask import Flask, jsonify, render_template_string, request


inTitle = ('Specificity Network Identification via '
           'Positional Entropy based Refinement (SNIPER)')


app = Flask(__name__)


@app.route('/upload', methods=['POST'])
def upload():
    fastaFile = request.files['fastaFile']
    # process the uploaded file...
    return f"Uploaded: {fastaFile.filename}"



@app.route('/run', methods=['POST'])
def run():
    try:
        # Get data from the request
        data = request.get_json()
        enzymeName = data.get('enzymeName')
        threshold = data.get('minS')
        NSubs = data.get('N')

        # Simulate some processing (replace this with your actual logic)
        result = {
            "message": "Processed successfully",
            "enzyme": enzymeName,
            "minS": threshold,
            "NSubs": NSubs
        }
        for key, value in result.items():
            print(key, value)

        # Return the result as a JSON response
        return jsonify(result)

    except Exception as e:
        # Return error message as JSON
        return jsonify({"error": f"Error: {str(e)}"}), 400



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
                    padding: {{ padTB }}px {{ padSide }}px {{ padTB }}px {{ padSide }}px;
                    margin: {{ spacer }}px auto;
                    margin-bottom: {{ marginB }}px;
                }
                .container2 {
                    background-color: {{ grey }};
                    font-size: 20px;
                    text-align: center;
                    border-radius: {{ borderRad }}px;
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    padding: {{ padTB }}px {{ padSide }}px {{ marginB }}px {{ padSide }}px;
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
                    margin-bottom: {{ marginB }}px;
                    border: 2px solid {{ greyDark }};
                    border-radius: {{ borderRad }}px;
                    width: 100%;
                    max-width: 350px;
                    padding: {{ padInput }}px;
                    
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
            <script>
                function evaluateForm() {
                    const enzymeName = document.getElementById('enzymeName').value;
                    const minS = document.getElementById('minS').value;
                    const N = document.getElementById('N').value;
                
                    fetch('/run', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json',
                        },
                        body: JSON.stringify({
                            enzymeName: enzymeName,
                            threshold: minS,
                            topN: N
                        })
                    })
                    .then(response => response.json())
                    .then(data => {
                        console.log("Response Data:", data);  // Log the data to ensure it's correct
                    
                        const newDiv = document.createElement('div');
                        newDiv.className = 'container';
                        const formattedJSON = JSON.stringify(data, null, 2); // pretty JSON
                    
                        newDiv.innerHTML = `
                            <div class="div-header">Results:</div>
                            <p><strong>Enzyme:</strong> ${enzymeName}</p>
                            <p><strong>Threshold:</strong> ${minS}</p>
                            <p><strong>Top N:</strong> ${N}</p>
                            <pre style="text-align: left; background-color: #1a1a1a; padding: 10px; border-radius: 8px; overflow-x: auto;">${formattedJSON}</pre>
                        `;
                        document.body.appendChild(newDiv);
                    })
                    .catch(error => {
                        console.error('Error:', error);
                        const errorDiv = document.createElement('div');
                        errorDiv.className = 'container';
                        errorDiv.innerHTML = `
                            <p style="color: red;">
                                Error fetching results: ${error.message || error}</p>`;
                        document.body.appendChild(errorDiv);
                    });

                }
                </script>
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
            <form class="container" method="POST" 
                    action="/upload" enctype="multipart/form-data">
                <div class="div-header">Upload FASTA File:</div>
                <input type="file" name="fastaFile" accept=".fasta,.fa" required
                class="input-form">
                <button type="submit">Upload</button>
            </form>
            <div class="container2">
                <form action="/run" method="POST" class="input-form">
                    <div class="div-header">Experiment Parameters:</div>
                    <div class="form-group">
                        <label for="enzymeName">Enzyme Name:</label>
                        <input type="text" id="enzymeName" name="enzymeName" required>
                    </div>

                    <div class="form-group">
                        <label for="minS">Entropy Threshold:</label>
                        <input type="number" id="minS" name="minS" 
                        step="0.1" required>
                    </div>

                    <div class="form-group">
                        <label for="N">Select Number of Substrates:</label>
                        <input type="number" id="N" name="N" required>
                    </div>
            
                    <button type="button" onclick="evaluateForm()">Evaluate</button>

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


    title=inTitle,

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
