from flask import Flask, jsonify, render_template, request
from functions import processData, subsDefault


app = Flask(__name__)


# Add citations
# Position descriptions for Trie
# Rename algorithm

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
            # Use template substrates
            substrates = subsDefault()

    except Exception as e:
        return jsonify({"error": f"Error: {str(e)}"}), 400


    # Get other data from the form
    enzymeName = request.form.get('enzymeName') if loadFile else 'Template Enzyme'
    entropyMin = request.form.get('entropyMin')
    NSelect = int(request.form.get('N'))

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
        pg1="For enzymes that "
            " This program will take substrates for a given enzyme and identify the "
            "Motif, of the recognition sequence within the larger protein sequence. "
            "The Motif is identified by the positions in the substrate that have "
            "Entropy scores (∆S) that exceed a entropyMin value.",
        pg2="∆S is evaluated at each position in the substrate sequence and is found by "
            "the difference between the Maximum Entropy (S<sub>Max</sub>) and the "
            "Shannon Entropy (S<sub>Shannon</sub>)",
        equation="∆S = S<sub>Max</sub> - S<sub>Shannon</sub> = log<sub>2</sub>(20) - "
                "∑(-prob<sub>AA</sub> * log<sub>2</sub>(prob<sub>AA</sub>))",
        pg3="Once we have identified the Motif, we can bin the substrates to select the "
            "only the recognition sequence. In other words, we remove the parts of the "
            "substrate that are not important for an Enzyme-Substrate interaction.",
        pg4="We can count the occurrences of each Motif, and plot this data in a Bar "
            "Graph, and a Word Cloud. This allows us to display the observed "
            "combinations of amino acids that fit into the enzymes active site.",
        pg5="We can further evaluate the Motif by feeding the data into a Suffix Tree. "
            "This analysis will select only the important residues within the motif, "
            "and plot the amino acids as nodes with lines connecting the observed "
            "combinations. This reveals the Specificity Network of your enzyme, as we "
            "can visualize unique preferences for a given amino acid when another is "
            "present in a preferred substrate.",
        pg6="The source code is available at "
            "<a href=\"https://github.com/Collinformatics/Webpage\" "
            "target=\"_blank\">Github/Collinformatics</a>."
    )


if __name__ == '__main__':
    app.run(debug=True)
