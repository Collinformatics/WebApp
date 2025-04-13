from flask import Flask, jsonify, render_template, request
from functions import processData, subsDefault


app = Flask(__name__)


@app.route('/run', methods=['POST'])
def run():
    # Get file from the form
    try:
        substrates = request.files.get('textFile')
        loadFile = False
        if substrates:
            loadFile = True

            # Read the contents of the file
            substrate = substrates.read().decode('utf-8')

            # Split the content by lines to get each substrate sequence
            substrates = substrate.splitlines()
        else:
            substrates = subsDefault()

    except Exception as e:
        return jsonify({"error": f"Error A: {str(e)}"}), 400


    # Get other data from the form
    enzymeName = request.form.get('enzymeName')
    entropyMin = request.form.get('entropyMin')
    NSelect = request.form.get('N')
    NSelect = int(NSelect)

    # Evaluate: Data
    dataset = processData(substrates, entropyMin, NSelect, enzymeName, loadFile)

    result = {
        "enzyme": enzymeName,
        "entropyMin": entropyMin,
        "NSelect": NSelect,
        "NBinSubs": dataset['NBinSubs'],
        "figProb": dataset['probability'],
        "figEntropy": dataset['entropy'],
        "figLogo": dataset['pLogo'],
        "barCounts": dataset['barCounts'],
        "barProb": dataset['barProb'],
        "figWords": dataset['words'],
        "figTrie": dataset['trie']
    }

    return jsonify(result)


@app.route('/')
def home():
    return render_template('home.html',
        black='#151515', white='FFFFFF',
        grey='#303030', greyDark='#202020',
        green='#23FF55', greenHighlight='#1AD747',
        spacer=20, spacerMini=5,
        padSide=50, padTB=30, padInput=8,
        marginB=12, marginButton=12,
        fontSize=16,
        borderRad=5,
        title="Specificity Network Identification via Positional Entropy based Refinement "
              "(SNIPER)",
        pg1="This program will take substrates for a given enzyme and identify the "
            "Motif, of the recognition sequence within the larger protein sequence. "
            "The Motif is identified by the positions in the substrate that have "
            "Entropy scores (∆S) that exceed a entropyMin value.",
        pg2="∆S is evaluated at each position in the substrate sequence and is found by "
            "the difference between the Maximum Entropy (S<sub>Max</sub>) and the "
            "Shannon Entropy (S<sub>Shannon</sub>)",
        equation="∆S = S<sub>Max</sub> - S<sub>Shannon</sub> = log<sub>2</sub>(20) - "
                "∑(-prob<sub>AA</sub> * log<sub>2</sub>(prob<sub>AA</sub>))",
        pg3="Once we have identified the Motif, we can bin the substrates to select the only"
            "the recognition sequence. In other words, we remove the parts of the substrate "
            "that are not important for an Enzyme-Substrate interaction.",
        pg4="We can count the occurrences of each Motif, and plot this data in a Bar Graph, "
            "and a Word Cloud. This allows us to display the optimal combinations of amino "
            "acids that fit into the enzymes active site."
            "will plot the and use them as input for a Suffix Tree. This will plot the "
            "preferred residues in descending order of Specificity, as determined by the "
            "∆S values.",
        pg5="We can further evaluate the Motif by feeding the data into a Suffix Tree. "
            "This will reveal the Specificity Network of your enzyme, or the specific "
            "connections between the preferences for a given amino acid when another is "
            "present."
    )



if __name__ == '__main__':
    app.run(debug=True)
