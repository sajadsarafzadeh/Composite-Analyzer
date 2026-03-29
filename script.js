let layers = [];

// ================= MATERIAL DATABASE =================
const materials = {
    carbon: {
        name: "Carbon/Epoxy",
        E1: 181e9,
        E2: 10.3e9,
        v12: 0.28,
        G12: 7.17e9
    },
    glass: {
        name: "Glass/Epoxy",
        E1: 38.6e9,
        E2: 8.27e9,
        v12: 0.26,
        G12: 4.14e9
    },
    kevlar: {
        name: "Kevlar/Epoxy",
        E1: 76e9,
        E2: 5.5e9,
        v12: 0.34,
        G12: 2.3e9
    }
};

// ================= UI =================

function addLayer() {
    const theta = parseFloat(document.getElementById('theta').value);
    const t = parseFloat(document.getElementById('thick').value) * 1e-3;
    const materialKey = document.getElementById('material').value;

    layers.push({
        theta,
        t,
        material: materials[materialKey]
    });

    updateTable();
}

function updateTable() {
    const tbody = document.querySelector('#layupTable tbody');

    tbody.innerHTML = layers.map((l, i) =>
        `<tr>
            <td>${i+1}</td>
            <td>${l.material.name}</td>
            <td>${l.theta}</td>
            <td>${l.t.toExponential(3)}</td>
        </tr>`
    ).join('');
}

// ================= CLT =================

function getQ(E1, E2, v12, G12) {
    const v21 = (v12 * E2) / E1;
    const denom = 1 - v12 * v21;

    return [
        [E1/denom, v12*E2/denom, 0],
        [v12*E2/denom, E2/denom, 0],
        [0, 0, G12]
    ];
}

function transformQ(Q, theta) {
    const rad = theta * Math.PI / 180;
    const m = Math.cos(rad);
    const n = Math.sin(rad);

    const Q11 = Q[0][0];
    const Q22 = Q[1][1];
    const Q12 = Q[0][1];
    const Q66 = Q[2][2];

    return [
        [
            Q11*m**4 + Q22*n**4 + 2*(Q12+2*Q66)*m**2*n**2,
            (Q11+Q22-4*Q66)*m**2*n**2 + Q12*(m**4+n**4),
            (Q11-Q12-2*Q66)*m**3*n - (Q22-Q12-2*Q66)*m*n**3
        ],
        [
            (Q11+Q22-4*Q66)*m**2*n**2 + Q12*(m**4+n**4),
            Q11*n**4 + Q22*m**4 + 2*(Q12+2*Q66)*m**2*n**2,
            (Q11-Q12-2*Q66)*m*n**3 - (Q22-Q12-2*Q66)*m**3*n
        ],
        [
            (Q11-Q12-2*Q66)*m**3*n - (Q22-Q12-2*Q66)*m*n**3,
            (Q11-Q12-2*Q66)*m*n**3 - (Q22-Q12-2*Q66)*m**3*n,
            (Q11+Q22-2*Q12-2*Q66)*m**2*n**2 + Q66*(m**4+n**4)
        ]
    ];
}

function calculateABD() {

    if (layers.length === 0) {
        alert("Add at least one layer!");
        return;
    }

    // total thickness
    let h = layers.reduce((sum, l) => sum + l.t, 0);

    // z coordinates
    let z = [-h / 2];
    for (let i = 0; i < layers.length; i++) {
        z.push(z[i] + layers[i].t);
    }

    let A = [[0,0,0],[0,0,0],[0,0,0]];
    let B = [[0,0,0],[0,0,0],[0,0,0]];
    let D = [[0,0,0],[0,0,0],[0,0,0]];

    layers.forEach((layer, k) => {

        const mat = layer.material;

        const Q = getQ(mat.E1, mat.E2, mat.v12, mat.G12);
        const Qbar = transformQ(Q, layer.theta);

        const z0 = z[k];
        const z1 = z[k+1];

        for (let i=0;i<3;i++){
            for (let j=0;j<3;j++){
                A[i][j] += Qbar[i][j]*(z1 - z0);
                B[i][j] += 0.5*Qbar[i][j]*(z1**2 - z0**2);
                D[i][j] += (1/3)*Qbar[i][j]*(z1**3 - z0**3);
            }
        }
    });

    displayMatrix("a-matrix", A);
    displayMatrix("b-matrix", B);
    displayMatrix("d-matrix", D);

    console.log("A:", A);
    console.log("B:", B);
    console.log("D:", D);
}

// ================= DISPLAY =================

function displayMatrix(id, M) {
    const el = document.getElementById(id);

    let html = "<table class='matrix'>";
    M.forEach(row => {
        html += "<tr>" + row.map(v => `<td>${v.toExponential(3)}</td>`).join("") + "</tr>";
    });
    html += "</table>";

    el.innerHTML = html;
}
