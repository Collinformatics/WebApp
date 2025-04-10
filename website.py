from flask import Flask, render_template_string

inTitle = ('Specificity Network Identification via<br>'
           'Positional Entropy based Refinement')

app = Flask(__name__)

@app.route('/')
def home():
    return render_template_string('''
    <html>
        <head>
            <title>SNIPER</title>
            <style>
                body {
                    padding: 40px;
                }
                h1 {
                    text-align: center;
                    font-size: 35px;
                    padding-bottom: 15px;
                }
                p {
                    font-size: 25px;
                    padding-top: 20px;
                    padding-left: 20px;
                }
                .custom-div {
                    margin-top: 30px;
                    padding: 20px;
                    background-color: #f2f2f2;
                    border: 1px solid #ccc;
                    font-size: 20px;
                }
            </style>
        </head>
        <body>
            <h1>{{ title|safe }}</h1>
            <p>{{ content }}</p>
            <div class="custom-div">
                New div
            </div>
        </body>
    </html>
    ''', title=inTitle, content="Words")

if __name__ == '__main__':
    app.run(debug=True)
