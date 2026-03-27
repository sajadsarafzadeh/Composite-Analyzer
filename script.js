let layers = [];

// Material Constants (Graphite/Epoxy example from your code)
const E1 = 181e9, E2 = 10.3e9, v12 = 0.28, G12 = 7.17e9;
const v21 = v12 / (E1 / E2);

function addLayer() {
    const theta = parseFloat(document.getElementById('theta').value);
    const t = parseFloat(document.getElementById('thick').value) * 1e-3; // mm to m
    layers.push({ theta, t });
    updateTable();
}

function updateTable() {
    const tbody = document.querySelector('#layupTable tbody');
    tbody.innerHTML = layers.map((l, i) => 
        `<tr><td>${i+1}</td><td>${l.theta}</td><td>${l.t.toFixed(5)}</td></tr>`
    ).join('');
}

function calculateABD() {
    let totalThick = layers.reduce((sum, l) => sum + l.t, 0);
    let z = [-totalThick / 2];
    let currentZ = -totalThick / 2;
    
    layers.forEach(l => {
        currentZ += l.t;
        z.push(currentZ);
    });

    let A = [[0,0,0],[0,0,0],[0,0,0]];
    
    // Q Matrix (On-axis)
    const Q11 = E1 / (1 - v12 * v21);
    const Q12 = (v21 * E1) / (1 - v12 * v21);
    const Q22 = E2 / (1 - v12 * v21);
    const Q66 = G12;

    layers.forEach((layer, i) => {
        const rad = layer.theta * Math.PI / 180;
        const m = Math.cos(rad);
        const n = Math.sin(rad);

        // Q-bar Transformation (Off-axis)
        const Qbar11 = Q11*Math.pow(m,4) + 2*(Q12+2*Q66)*Math.pow(m,2)*Math.pow(n,2) + Q22*Math.pow(n,4);
        // ... (Other Qbar components would be added here similar to your MATLAB U-parameters)

        // Simple A-matrix sum: A = sum(Qbar * thickness)
        A[0][0] += Qbar11 * layer.t;
    });

    document.getElementById('a-matrix').innerText = JSON.stringify(A[0][0].toExponential(3));
    alert("A[1,1] Calculated! (Full matrix logic follows same loop)");
}
