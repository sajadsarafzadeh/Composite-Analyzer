let layers = [];

const materials = {
    carbon: { 
        name: "Carbon/Epoxy", 
        E1: 181e9, E2: 10.3e9, v12: 0.28, G12: 7.17e9,
        F_xx: 0.444e-18, F_yy: 101.6e-18, F_xy: -3.36e-18, F_ss: 216.2e-18, F_x: 0, F_y: 20.93e-9
    },
    glass: { 
        name: "Glass/Epoxy", 
        E1: 38.6e9, E2: 8.27e9, v12: 0.26, G12: 4.14e9,
        F_xx: 0.25e-18, F_yy: 50.0e-18, F_xy: -2.0e-18, F_ss: 120.0e-18, F_x: 0, F_y: 10.0e-9
    },
    kevlar: { 
        name: "Kevlar/Epoxy", 
        E1: 76e9, E2: 5.5e9, v12: 0.34, G12: 2.3e9,
        F_xx: 0.5e-18, F_yy: 80.0e-18, F_xy: -3.0e-18, F_ss: 200.0e-18, F_x: 0, F_y: 15.0e-9
    }
};

function addLayer() {
    const theta = parseFloat(document.getElementById('theta').value);
    const t = parseFloat(document.getElementById('thick').value) * 1e-3;
    const mat = materials[document.getElementById('material').value];
    layers.push({theta, t, material: mat});
    updateTable();
}

function updateTable() {
    const tbody = document.querySelector('#layupTable tbody');
    tbody.innerHTML = layers.map((l,i)=>`
        <tr>
            <td>${i+1}</td>
            <td>${l.material.name}</td>
            <td>${l.theta}</td>
            <td>${l.t.toExponential(2)}</td>
        </tr>`).join('');
}

function getQ(m) {
    const v21 = (m.v12 * m.E2)/m.E1;
    const d = 1 - m.v12*v21;
    return [
        [m.E1/d, m.v12*m.E2/d, 0],
        [m.v12*m.E2/d, m.E2/d, 0],
        [0, 0, m.G12]
    ];
}

function Qbar(Q, theta) {
    const m = Math.cos(theta*Math.PI/180);
    const n = Math.sin(theta*Math.PI/180);
    const Q11=Q[0][0], Q22=Q[1][1], Q12=Q[0][1], Q66=Q[2][2];
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
    if(layers.length === 0) { alert("Add layers first!"); return; }

    // Thickness & z coordinates
    let h = layers.reduce((sum,l)=>sum + l.t, 0);
    let z = [-h/2];
    layers.forEach(layer => { z.push(z[z.length-1] + layer.t); });

    // Initialize matrices
    let A = math.zeros(3,3)._data;
    let B = math.zeros(3,3)._data;
    let D = math.zeros(3,3)._data;

    // Loop layers
    layers.forEach((layer,k)=>{
        const Q = getQ(layer.material);
        const Qb = Qbar(Q, layer.theta);
        const z0 = z[k], z1 = z[k+1];
        for(let i=0;i<3;i++){
            for(let j=0;j<3;j++){
                A[i][j] += Qb[i][j]*(z1-z0);
                B[i][j] += 0.5*Qb[i][j]*(z1**2 - z0**2);
                D[i][j] += (1/3)*Qb[i][j]*(z1**3 - z0**3);
            }
        }
    });

    // Display A B D
    displayMatrix('A',A);
    displayMatrix('B',B);
    displayMatrix('D',D);

    // Build ABD and inverse
    const ABD = [
        [...A[0], ...B[0]],
        [...A[1], ...B[1]],
        [...A[2], ...B[2]],
        [...B[0], ...D[0]],
        [...B[1], ...D[1]],
        [...B[2], ...D[2]]
    ];
    const ABD_inv = math.inv(ABD);

    displayMatrix("ABD", ABD);
    displayMatrix("ABD_inv", ABD_inv);

    // Get loads
    const NM = [
        +document.getElementById("N1").value,
        +document.getElementById("N2").value,
        +document.getElementById("N6").value,
        +document.getElementById("M1").value,
        +document.getElementById("M2").value,
        +document.getElementById("M6").value
    ];

    // Solve ε0 & κ
    const result = math.multiply(ABD_inv, NM);
    const eps0 = result.slice(0,3);
    const kappa = result.slice(3,6);

    // Strain at z helper
    function strainAtZ(zval){ return [eps0[0]+kappa[0]*zval, eps0[1]+kappa[1]*zval, eps0[2]+kappa[2]*zval]; }

    // Stress & Tsai-Wu FI
    let resultHTML = "";
    layers.forEach((layer,k)=>{
        const Q = getQ(layer.material);
        const Qb = Qbar(Q, layer.theta);
        const zmid = (z[k]+z[k+1])/2;
        const strain = strainAtZ(zmid);
        const stress = math.multiply(Qb, strain);
        const m = layer.material;
        let FI = 0;
        if(m.F_xx){
            FI = m.F_xx*stress[0]**2 + m.F_yy*stress[1]**2 + m.F_ss*stress[2]**2 + 2*m.F_xy*stress[0]*stress[1] + m.F_x*stress[0] + m.F_y*stress[1];
        }
        resultHTML += `
            Layer ${k+1} (${m.name}) <br>
            σ = [${stress[0].toExponential(2)}, ${stress[1].toExponential(2)}, ${stress[2].toExponential(2)}] <br>
            FI = ${FI.toFixed(3)} → ${FI>=1 ? "FAIL ❌" : "SAFE ✅"} 
            <hr>
        `;
    });
    document.getElementById("failure").innerHTML = resultHTML;
}

function displayMatrix(id,M){
    document.getElementById(id).innerHTML =
        "<table>"+M.map(r=>"<tr>"+r.map(v=>`<td>${v.toExponential(2)}</td>`).join("")+"</tr>").join("")+"</table>";
}
