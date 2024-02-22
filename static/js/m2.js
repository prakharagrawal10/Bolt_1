document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('imageForm').addEventListener('submit', function(event) {
        event.preventDefault(); // Prevent default form submission
        
        var input = document.getElementById('imageFile');
        var file = input.files[0];

        var reader = new FileReader();
        reader.onload = function(event) {
            var imageData = event.target.result;

            // Make AJAX request to Flask server
            var xhr = new XMLHttpRequest();
            xhr.open('POST', '/predict_model2');
            xhr.setRequestHeader('Content-Type', 'application/json');
            xhr.onload = function() {
                if (xhr.status === 200) {
                    var result = JSON.parse(xhr.responseText);
                    displayVisualization(result);
                } else {
                    console.error('Error:', xhr.statusText);
                }
            };
            xhr.onerror = function() {
                console.error('Request failed');
            };
            
            // Send the image data as base64 string
            xhr.send(JSON.stringify({ imageData: imageData }));
        };

        // Read the selected file as data URL (base64)
        reader.readAsDataURL(file);
    });

    function displayVisualization(result) {
        // Display the visualization using the result data
        // Implement this according to your requirements
    }
});
