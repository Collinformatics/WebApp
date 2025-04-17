from flask import Flask, jsonify, render_template, request
from functions import Webpage


app = Flask(__name__)

# Initialize: Application
web = Webpage()


@app.route('/run', methods=['POST'])
def run():
    try:
        message = request.files.get('message')
        print(f'Received message: {message}')

        web.buttonState()

    except Exception as e:
        return jsonify({"error": f"Error: {str(e)}"}), 400




    # result = {
    #     "message": enzymeName,
    #     "entropyMin": entropyMin,
    # }
    #
    # return jsonify(result)


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
    )


if __name__ == '__main__':
    app.run(debug=True)
