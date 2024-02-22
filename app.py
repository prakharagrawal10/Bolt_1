from flask import Flask, render_template, request, jsonify
from model import get_molecular_properties
from model2 import process_image_data  # Import function from model2.py

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('services.html')

@app.route('/m1')
def model1():
    return render_template('m1.html')

@app.route('/m2')
def model2():
    return render_template('m2.html')

@app.route('/m3')
def model3():
    return render_template('Cancer_data_analsis.html')

@app.route('/predict', methods=['POST'])
def predict():
    data = request.json
    if data is None or 'inputData' not in data:
        return jsonify({'error': 'Invalid request. Missing input data.'}), 400

    smiles = data['inputData']
    if not smiles:
        return jsonify({'error': 'Empty input data.'}), 400

    try:
        result = get_molecular_properties(smiles)
        return jsonify(result), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/predict_model2', methods=['POST'])
def predict_model2():
    print("Hello world! Submit button pressed.")  # Add this line

    data = request.json
    if data is None or 'imageData' not in data:
        return jsonify({'error': 'Invalid request. Missing input data.'}), 400

    image_data = data['imageData']
    if not image_data:
        return jsonify({'error': 'Empty input data.'}), 400

    try:
        result = process_image_data(image_data)
        return jsonify(result), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)

