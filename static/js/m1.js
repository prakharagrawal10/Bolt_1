document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('predictionForm').addEventListener('submit', function(event) {
        event.preventDefault(); // Prevent default form submission

        var smiles = document.getElementById('inputSmiles').value;

        // Make an AJAX request to Flask server
        var xhr = new XMLHttpRequest();
        xhr.open('POST', '/predict');
        xhr.setRequestHeader('Content-Type', 'application/json');
        xhr.onload = function() {
            if (xhr.status === 200) {
                var result = JSON.parse(xhr.responseText);
                displayResult(result);
            } else {
                console.error('Error:', xhr.statusText);
            }
        };
        xhr.onerror = function() {
            console.error('Request failed');
        };
        xhr.send(JSON.stringify({ inputData: smiles }));
    });

    function displayResult(result) {
        var resultDiv = document.getElementById('result');
        resultDiv.innerHTML = ''; // Clear previous result
        if (result.error) {
            resultDiv.innerText = 'Error: ' + result.error;
        } else {
            // Display molecular properties
            resultDiv.innerHTML = `
                <p>Number of atoms: ${result.num_atoms}</p>
                <p>Number of bonds: ${result.num_bonds}</p>
                <p>Molecular weight: ${result.molecular_weight.toFixed(2)} g/mol</p>
                <p>CLogP (lipophilicity): ${result.logP.toFixed(2)}</p>
                <p>Number of rotatable bonds: ${result.num_rotatable_bonds}</p>
                <p>Hydrogen bond donors: ${result.hydrogen_bond_donors}</p>
                <p>Hydrogen bond acceptors: ${result.hydrogen_bond_acceptors}</p>
                <p>Polar surface area: ${result.polar_surface_area.toFixed(2)} Å^2</p>
                <p>Molecular formula: ${result.molecular_formula}</p>
                <p>Aromatic atoms: ${result.aromatic_atoms}</p>
                <p>Ring count: ${result.ring_count}</p>
                <img src="${result.img_path}" alt="Molecule">
            `;
        }
    }
});
