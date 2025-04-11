from flask import Flask, render_template_string

# <br>
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
                    background-color: #151515;
                    color: #f0f0f0;
                    font-family: Arial, sans-serif;
                    padding: 40px;
                }
                h1 {
                    text-align: center;
                    font-size: 35px;
                    padding-bottom: 10px;
                }
                p {
                    font-size: 20px;
                    text-align: center;
                    padding-top: 5px;
                    padding-bottom: 0px;
                }
                .upload-container {
                    background-color: #303030;
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    border-radius: 10px;
                    width: fit-content;
                    margin: 0px auto;
                    box-shadow: 0 0 10px rgba(0, 0, 0, 0.5);
                    padding: 30px 50px 30px 50px;
                }
                .upload-header {
                    color: #f0f0f0;
                    font-size: 25px;
                    font-weight: bold;
                    padding-top: 5px;
                    padding-bottom: 20px;
                }
                input[type="file"] {
                    background-color: #151515;
                    color: #101010;
                    font-size: 18px;
                    margin-bottom: 15px;
                    border-radius: 6px;
                    width: 100%;
                    max-width: 350px;
                    border: 2px solid #202020;
                    padding: 8px;

                }
                button {
                    background-color: #24E100;
                    color: #FFFFFF;
                    font-size: 18px;
                    cursor: pointer;
                    border: none;
                    border-radius: 5px;
                    padding: 10px 20px;
                }
                button:hover {
                    background-color: #17FF44;
                }
                .custom-div {
                    background-color: #303030;
                    color: #EE1088;
                    font-size: 20px;
                    margin-top: 30px;
                    border-radius: 5px;
                    text-align: center;
                    padding: 20px;
                }
            </style>
        </head>
        <body>
            <h1>{{ title|safe }}</h1>
            <p>{{ ct1|safe }}</p>
            <p>{{ ct2|safe }}</p>
            <p>{{ equation|safe }}</p>
            <p>{{ ct3|safe }}</p>
            <form class="upload-container" method="POST" 
                    action="/upload" enctype="multipart/form-data">
                <div class="upload-header">Upload FASTA File:</div>
                <input type="file" name="fastaFile" accept=".fasta,.fa" required>
                <button type="submit">Upload</button>
            </form>
            <div class="custom-div">
                New div
            </div>
        </body>
    </html>
    ''',
    title=inTitle,

    ct1="This program will take substrates for a given enzyme and identify the "
        "Motif, of the recognition sequence within the larger protein sequence. "
        "The Motif is identified by the positions in the substrate that have "
        "Entropy scores (∆S) that exceed a threshold value.",

    ct2="∆S is evaluated at each position in the substrate sequence and is found by "
        "the difference between the Maximum Entropy (S<sub>Max</sub>) and the "
        "Shannon Entropy (S<sub>Shannon</sub>)",

    equation="∆S = S<sub>Max</sub> - S<sub>Shannon</sub> = log<sub>2</sub>(20) - "
            "∑(-prob<sub>AA</sub> * log<sub>2</sub>(prob<sub>AA</sub>))",

    ct3="")



if __name__ == '__main__':
    app.run(debug=True)
