from flask import Flask, jsonify, render_template, request
from functions import WebApp


app = Flask(__name__)

# Initialize: Application
web = WebApp()


@app.route('/run', methods=['POST'])
def run():
    try:
        message = request.form.get('message')
        messages = web.getMessage(message)

    except Exception as e:
        return jsonify({"error": f"Error: {str(e)}"}), 400

    result = {
        "message": message,
        "messages": messages,
    }

    return jsonify(result)


@app.route('/')
def home():
    return render_template('home.html',
        black='#151515', white='#FFFFFF',
        colorBG='#303030', colorAccent='#454545',
        colorHeader='#23FF55', buttonHighlight='#1AD747',
        spacer=20, spacerMini=5,
        padSide=25, padTB=30, padInput=8,
        marginB=12, marginButton=12,
        fontSize=16,
        borderRad=5,
    )


if __name__ == '__main__':
    app.run(debug=True)
