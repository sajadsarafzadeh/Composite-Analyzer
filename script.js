let layers = [];

const materials = {
    carbon: { name:"Carbon/Epoxy", E1:181e9, E2:10.3e9, v12:0.28, G12:7.17e9 },
    glass: { name:"Glass/Epoxy", E1:38.6e9, E2:8.27e9, v12:0.26, G12:4.14e9 },
    kevlar: { name:"Kevlar/Epoxy", E1:76e9, E2:5.5e9, v12:0.34, G12:2.3e9 }
};

function addLayer(){
    const theta = parseFloat(document.getElementById('theta').value);
    const t = parseFloat(document.getElementById('thick').value) * 1e-3;
    const mat = materials[document.getElementById('material').value];

    layers.push({theta,t,material:mat});
    updateTable();
}

function updateTable(){
    const tbody = document.querySelector('#layupTable tbody');
    tbody.innerHTML = layers.map((l,i)=>`
        <tr>
            <td>${i+1}</td>
            <td>${l.material.name}</td>
            <td>${l.theta}</td>
            <td>${l.t.toExponential(2)}</td>
        </tr>`).join('');
}

function getQ(m){
    const v21 = (m.v12 * m.E2)/m.E1;
    const d = 1 - m.v12*v21;
    return [
        [m.E1/d, m.v12*m.E2/d, 0],
        [m.v12*m.E2/d, m.E2/d, 0],
        [0,0,m.G12]
    ];
}

function Qbar(Q, theta){
    const m = Math.cos(theta*Math.PI/180);
    const n = Math.sin(theta*Math.PI/180);
    const Q11=Q[0][0],Q22=Q[1][1],Q12=Q[0][1],Q66=Q[2][2];

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

function calculateABD(){

    if(layers.length === 0){
        alert("Add layers first!");
        return;
    }

    // ===== 1. Thickness & z =====
    let h = layers.reduce((sum,l)=>sum + l.t, 0);

    let z = [-h/2];
    layers.forEach(layer=>{
        z.push(z[z.length-1] + layer.t);
    });

    // ===== 2. Initialize Matrices =====
    let A = [[0,0,0],[0,0,0],[0,0,0]];
    let B = [[0,0,0],[0,0,0],[0,0,0]];
    let D = [[0,0,0],[0,0,0],[0,0,0]];

    // ===== 3. Loop Layers =====
    layers.forEach((layer, k)=>{

        const Q = getQ(layer.material);
        const Qb = Qbar(Q, layer.theta);

        let z0 = z[k];
        let z1 = z[k+1];

        for(let i=0;i<3;i++){
            for(let j=0;j<3;j++){

                A[i][j] += Qb[i][j] * (z1 - z0);
                B[i][j] += 0.5 * Qb[i][j] * (z1**2 - z0**2);
                D[i][j] += (1/3) * Qb[i][j] * (z1**3 - z0**3);

            }
        }
    });

    // ===== 4. Display A B D =====
    displayMatrix('A', A);
    displayMatrix('B', B);
    displayMatrix('D', D);

    // ===== 5. Build ABD =====
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
    displayMatrix("abbd", ABD_inv);

    // ===== 6. Loads =====
    const NM = [
        +document.getElementById("N1").value,
        +document.getElementById("N2").value,
        +document.getElementById("N6").value,
        +document.getElementById("M1").value,
        +document.getElementById("M2").value,
        +document.getElementById("M6").value
    ];

    // ===== 7. Solve ε0 & κ =====
    const result = math.multiply(ABD_inv, NM);

    const eps0 = result.slice(0,3);
    const kappa = result.slice(3,6);

    console.log("eps0:", eps0);
    console.log("kappa:", kappa);

    // ===== 8. Helper =====
    function strainAtZ(zval){
        return [
            eps0[0] + kappa[0]*zval,
            eps0[1] + kappa[1]*zval,
            eps0[2] + kappa[2]*zval
        ];
    }

    // ===== 9. Stress + Failure =====
    let resultHTML = "";

    layers.forEach((layer, k)=>{

        const Q = getQ(layer.material);
        const Qb = Qbar(Q, layer.theta);

        const zmid = (z[k] + z[k+1]) / 2;

        const strain = strainAtZ(zmid);
        const stress = math.multiply(Qb, strain);

        // Tsai-Wu
        const s1 = stress[0];
        const s2 = stress[1];
        const s12 = stress[2];

        const m = layer.material;

        let FI = 0;

        if(m.F_xx){ // اگر تعریف شده بود
            FI =
                m.F_xx*s1*s1 +
                m.F_yy*s2*s2 +
                m.F_ss*s12*s12 +
                2*m.F_xy*s1*s2 +
                m.F_x*s1 +
                m.F_y*s2;
        }

        resultHTML += `
            Layer ${k+1} (${layer.material.name}) <br>
            σ = [${s1.toExponential(2)}, ${s2.toExponential(2)}, ${s12.toExponential(2)}] <br>
            FI = ${FI.toFixed(3)} → ${FI>=1 ? "FAIL ❌" : "SAFE ✅"} 
            <hr>
        `;
    });

    document.getElementById("failure").innerHTML = resultHTML;

    // ===== Debug =====
    console.log("A:", A);
    console.log("B:", B);
    console.log("D:", D);
}


function displayMatrix(id,M){
    document.getElementById(id).innerHTML =
        "<table>"+M.map(r=>"<tr>"+r.map(v=>`<td>${v.toExponential(2)}</td>`).join("")+"</tr>").join("")+"</table>";
}
