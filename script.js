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
        [Q11*m**4+2*(Q12+2*Q66)*m**2*n**2+Q22*n**4,0,0],
        [0,Q11*n**4+2*(Q12+2*Q66)*m**2*n**2+Q22*m**4,0],
        [0,0,Q66*(m**4+n**4)+(Q11+Q22-2*Q12-2*Q66)*m**2*n**2]
    ];
}

function calculateA(){
    let A=[[0,0,0],[0,0,0],[0,0,0]];
    layers.forEach(l=>{
        const Qb = Qbar(getQ(l.material), l.theta);
        for(let i=0;i<3;i++){
            for(let j=0;j<3;j++){
                A[i][j] += Qb[i][j]*l.t;
            }
        }
    });
    displayMatrix('A',A);
}

function displayMatrix(id,M){
    document.getElementById(id).innerHTML =
        "<table>"+M.map(r=>"<tr>"+r.map(v=>`<td>${v.toExponential(2)}</td>`).join("")+"</tr>").join("")+"</table>";
}
