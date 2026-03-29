let layers = [];

// ================= MATERIAL DATABASE =================
const materials = {
    carbon: {
        name: "Carbon/Epoxy",
        E1: 181e9, E2: 10.3e9, v12: 0.28, G12: 7.17e9,
        F_xx: 0.444e-18, F_yy: 101.6e-18, F_xy: -3.36e-18,
        F_ss: 216.2e-18, F_x: 0, F_y: 20.93e-9
    },
    glass: {
        name: "Glass/Epoxy",
        E1: 38.6e9, E2: 8.27e9, v12: 0.26, G12: 4.14e9,
        F_xx: 1.543e-18, F_yy: 273.3e-18, F_xy: -10.27e-18,
        F_ss: 192.9e-18, F_x: -0.697e-9, F_y: 23.78e-9
    },
    kevlar: {
        name: "Kevlar/Epoxy",
        E1: 76e9, E2: 5.5e9, v12: 0.34, G12: 2.3e9,
        F_xx: 3.039e-18, F_yy: 1572e-18, F_xy: -34.56e-18,
        F_ss: 865e-18, F_x: -3.541e-9, F_y: 64.46e-9
    }
};

// ================= UI =================

function addLayer(){
    const theta = parseFloat(theta.value);
    const t = parseFloat(thick.value) * 1e-3;
    const mat = materials[material.value];

    layers.push({theta, t, material: mat});
    updateTable();
}

function updateTable(){
    layupTable.querySelector("tbody").innerHTML =
        layers.map((l,i)=>`
        <tr>
            <td>${i+1}</td>
            <td>${l.material.name}</td>
            <td>${l.theta}</td>
            <td>${l.t.toExponential(2)}</td>
        </tr>`).join('');
}

// ================= CLT =================

function getQ(m){
    const v21 = (m.v12 * m.E2) / m.E1;
    const d = 1 - m.v12 * v21;

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
            Q11*m**4+Q22*n**4+2*(Q12+2*Q66)*m**2*n**2,
            (Q11+Q22-4*Q66)*m**2*n**2+Q12*(m**4+n**4),
            (Q11-Q12-2*Q66)*m**3*n-(Q22-Q12-2*Q66)*m*n**3
        ],
        [
            (Q11+Q22-4*Q66)*m**2*n**2+Q12*(m**4+n**4),
            Q11*n**4+Q22*m**4+2*(Q12+2*Q66)*m**2*n**2,
            (Q11-Q12-2*Q66)*m*n**3-(Q22-Q12-2*Q66)*m**3*n
        ],
        [
            (Q11-Q12-2*Q66)*m**3*n-(Q22-Q12-2*Q66)*m*n**3,
            (Q11-Q12-2*Q66)*m*n**3-(Q22-Q12-2*Q66)*m**3*n,
            (Q11+Q22-2*Q12-2*Q66)*m**2*n**2+Q66*(m**4+n**4)
        ]
    ];
}

// ================= MAIN =================

function analyze(){

    let h = layers.reduce((s,l)=>s+l.t,0);
    let z=[-h/2];
    layers.forEach(l=>z.push(z[z.length-1]+l.t));

    let A=[[0,0,0],[0,0,0],[0,0,0]];
    let B=[[0,0,0],[0,0,0],[0,0,0]];
    let D=[[0,0,0],[0,0,0],[0,0,0]];

    layers.forEach((l,k)=>{
        let Qb = Qbar(getQ(l.material), l.theta);
        let z0=z[k], z1=z[k+1];

        for(let i=0;i<3;i++){
            for(let j=0;j<3;j++){
                A[i][j]+=Qb[i][j]*(z1-z0);
                B[i][j]+=0.5*Qb[i][j]*(z1**2-z0**2);
                D[i][j]+=(1/3)*Qb[i][j]*(z1**3-z0**3);
            }
        }
    });

    display("A",A);
    display("B",B);
    display("D",D);

    // ===== LOADS =====
    const NM=[
        +N1.value,+N2.value,+N6.value,
        +M1.value,+M2.value,+M6.value
    ];

    const ABD=[
        [...A[0],...B[0]],
        [...A[1],...B[1]],
        [...A[2],...B[2]],
        [...B[0],...D[0]],
        [...B[1],...D[1]],
        [...B[2],...D[2]]
    ];

    const x = math.multiply(math.inv(ABD), NM);
    const eps0=x.slice(0,3);
    const k=x.slice(3,6);

    // ===== FAILURE =====
    let out="";

    layers.forEach((l,kp)=>{
        let zmid=(z[kp]+z[kp+1])/2;

        let strain=[
            eps0[0]+k[0]*zmid,
            eps0[1]+k[1]*zmid,
            eps0[2]+k[2]*zmid
        ];

        let stress=math.multiply(Qbar(getQ(l.material),l.theta),strain);

        let F=l.material;

        let FI =
            F.F_xx*stress[0]**2 +
            F.F_yy*stress[1]**2 +
            F.F_ss*stress[2]**2 +
            2*F.F_xy*stress[0]*stress[1] +
            F.F_x*stress[0] +
            F.F_y*stress[1];

        out += `Layer ${kp+1}: ${FI.toFixed(3)} → ${FI>=1?"FAIL ❌":"SAFE ✅"}<br>`;
    });

    failure.innerHTML=out;
}

// ================= DISPLAY =================

function display(id,M){
    document.getElementById(id).innerHTML =
        "<table>"+M.map(r=>"<tr>"+r.map(v=>`<td>${v.toExponential(2)}</td>`).join("")+"</tr>").join("")+"</table>";
}
