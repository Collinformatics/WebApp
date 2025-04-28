class WebApp:
    def __init__(self):
        self.buttonState = False
        self.messages = []

    def pressButton(self, message):
        print(f'Received data: {message}')

        return {'key': 'Returned data'}
