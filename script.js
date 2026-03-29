let layers = [];

const materials = {
    carbon: { 
        name: "Carbon/Epoxy", E1: 181e9, E2: 10.3e9, v12: 0.28, G12: 7.17e9,
        F_xx: 0.444e-18, F_yy: 101.6e-18, F_xy: -3.36e-18, F_ss: 216.2e-18, F_x:0, F_y:20.93e-9
    },
    glass: { 
        name: "Glass/Epoxy", E1: 38.6e9, E2: 8.27e9, v12: 0.26, G12: 4.14e9,
        F_xx:0.25e-18, F_yy:50e-18, F_xy:-2e-18, F_ss:120e-18, F_x:0, F_y:10e-9
    },
    kevlar: { 
        name: "Kevlar/Epoxy", E1: 76e9, E2: 5.5e9, v12: 0.34, G12: 2.3e9,
        F_xx:0.5e-18, F_yy:80e-18, F_xy:-3e-18, F_ss:200e-18, F_x:0, F_y:15e-9
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
        [Q11*m**4 + Q22*n**4 + 2*(Q12+2*Q66)*m**2*n**2,
         (Q11+Q22-4*Q66)*m**2*n**2 + Q12*(m**4+n**4),
         (Q11-Q12-2*Q66)*m**3*n - (Q22-Q12-2*Q66)*m*n**3],
        [(Q11+Q22-4*Q66)*m**2*n**2 + Q12*(m**4+n**4),
         Q11*n**4 + Q22*m**4 + 2*(Q12+2*Q66)*m**2*n**2,
         (Q11-Q12-2*Q66)*m*n**3 - (Q22-Q12-2*Q66)*m**3*n],
        [(Q11-Q12-2*Q66)*m**3*n - (Q22-Q12-2*Q66)*m*n**3,
         (Q11-Q12-2*Q66)*m*n**3 - (Q22-Q12-2*Q66)*m**3*n,
         (Q11+Q22-2*Q12-2*Q66)*m**2*n**2 + Q66*(m**4+n**4)]
    ];
}

function calculateABD() {
    if(layers.length===0){ alert("Add layers first!"); return; }

    let h = layers.reduce((sum,l)=>sum + l.t,0);
    let z = [-h/2];
    layers.forEach(l=>z.push(z[z.length-1]+l.t));

    let A = math.zeros(3,3)._data;
    let B = math.zeros(3,3)._data;
    let D = math.zeros(3,3)._data;

    layers.forEach((layer,k)=>{
        const Q = getQ(layer.material);
        const Qb = Qbar(Q, layer.theta);
        let z0=z[k], z1=z[k+1];
        for(let i=0;i<3;i++){for(let j=0;j<3;j++){
            A[i][j]+=Qb[i][j]*(z1-z0);
            B[i][j]+=0.5*Qb[i][j]*(z1**2-z0**2);
            D[i][j]+=(1/3)*Qb[i][j]*(z1**3-z0**3);
        }}
    });

    displayMatrix('A',A);
    displayMatrix('B',B);
    displayMatrix('D',D);

    const ABD = [
        [...A[0],...B[0]],
        [...A[1],...B[1]],
        [...A[2],...B[2]],
        [...B[0],...D[0]],
        [...B[1],...D[1]],
        [...B[2],...D[2]]
    ];

    const ABD_inv = math.inv(ABD);
    displayMatrix('ABD',ABD);
    displayMatrix('abbd',ABD_inv);

    // --- Tsai-Wu ---
    let eps0 = [0,0,0]; // initial, will be solved if loads provided
    let kappa = [0,0,0];
    let z_mid=[], strain_local=[], stress_local=[];
    let strain_global=[], stress_global=[], tsai=[];

    layers.forEach((layer,k)=>{
        const Q = getQ(layer.material);
        const Qb = Qbar(Q, layer.theta);
        const zmid=(z[k]+z[k+1])/2;
        z_mid.push(zmid);

        const eps_local=[...eps0]; // No external load, eps0=0
        const sigma_local=math.multiply(Qb,eps_local);

        const m = layer.material;
        const FI=m.F_xx*sigma_local[0]**2 + m.F_yy*sigma_local[1]**2 + m.F_ss*sigma_local[2]**2 + 2*m.F_xy*sigma_local[0]*sigma_local[1] + m.F_x*sigma_local[0] + m.F_y*sigma_local[1];
        tsai.push(FI);

        strain_local.push(eps_local);
        stress_local.push(sigma_local);
        // For global, rotate back by theta
        const c=Math.cos(layer.theta*Math.PI/180), s=Math.sin(layer.theta*Math.PI/180);
        const T = [[c*c,s*s,2*c*s],[s*s,c*c,-2*c*s],[-c*s,c*s,c*c-s*s]]; // simple approx
        strain_global.push(math.multiply(T,eps_local));
        stress_global.push(math.multiply(T,sigma_local));
    });

    // --- Plot ---
    let trace_eps1={x:z_mid,y:strain_global.map(e=>e[0]),name:"ε₁ (Global)",mode:"lines+markers"};
    let trace_eps2={x:z_mid,y:strain_global.map(e=>e[1]),name:"ε₂ (Global)",mode:"lines+markers"};
    let trace_eps12={x:z_mid,y:strain_global.map(e=>e[2]),name:"γ₁₂ (Global)",mode:"lines+markers"};

    let trace_sig1={x:z_mid,y:stress_global.map(s=>s[0]),name:"σ₁ (Global)",mode:"lines+markers"};
    let trace_sig2={x:z_mid,y:stress_global.map(s=>s[1]),name:"σ₂ (Global)",mode:"lines+markers"};
    let trace_sig12={x:z_mid,y:stress_global.map(s=>s[2]),name:"τ₁₂ (Global)",mode:"lines+markers"};

    Plotly.newPlot('plots',[trace_eps1,trace_eps2,trace_eps12,trace_sig1,trace_sig2,trace_sig12],
        {title:"Stress & Strain Profiles", xaxis:{title:"z (m)"}, yaxis:{title:"Value"}});
}

function displayMatrix(id,M){
    document.getElementById(id).innerHTML="<table>"+M.map(r=>"<tr>"+r.map(v=>`<td>${v.toExponential(2)}</td>`).join("")+"</tr>").join("")+"</table>";
}
