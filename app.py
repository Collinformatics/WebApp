from flask import Flask, jsonify, render_template, request
from functions import WebApp


app = Flask(__name__)

# Initialize: Application
webapp = WebApp()


@app.route('/run', methods=['POST'])
def run():
    data = request.get_json() # Get the JSON body
    message = data.get('message') # Extract the message

    # Call method
    data = webapp.pressButton(message)

    return jsonify(data)


@app.route('/')
def home():
    return render_template('home.html')


if __name__ == '__main__':
    app.run(debug=True)
