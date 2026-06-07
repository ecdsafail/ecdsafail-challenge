#!/usr/bin/env python3
"""Interactive qubit-allocation explorer: toggle stacked-by-role <-> scratch-only,
minimap window-slider, wheel/trackpad zoom (cursor-centered), drag-to-pan, peak markers."""
import json, sys
SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/qmap_frontier.tsv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/qmap_explorer.html"
TITLE = sys.argv[3] if len(sys.argv) > 3 else "frontier (1302q)"
rows, phases, pidx = [], [], {}
for line in open(SRC):
    if line.startswith("#") or not line.strip():
        continue
    c = line.rstrip("\n").split("\t")
    op, ph, active = int(c[0]), c[1], int(c[2])
    scr = int(c[3]) + int(c[4]); tx = int(c[5]) + int(c[6]); ty = int(c[7]) + int(c[8])
    u = int(c[9]) + int(c[10]); tr = int(c[11]) + int(c[12])
    if ph not in pidx:
        pidx[ph] = len(phases); phases.append(ph)
    rows.append([op, pidx[ph], active, scr, tx, ty, u, tr])
peak = max(r[2] for r in rows); maxscr = max(r[3] for r in rows)
D = json.dumps({"rows": rows, "phases": phases, "peak": peak, "maxscr": maxscr})
PAGE = r"""<!doctype html><html><head><meta charset=utf8><title>qmap explorer</title>
<style>body{font:13px -apple-system,system-ui,sans-serif;margin:0;background:#0f1115;color:#dde}
#wrap{padding:14px}h1{font-size:16px;margin:0 0 4px}
.bar{display:flex;gap:8px;align-items:center;margin:6px 0;flex-wrap:wrap}
button{background:#1a1e27;border:1px solid #333a48;color:#dde;padding:5px 11px;border-radius:6px;cursor:pointer}
button.on{background:#2a4d8f;border-color:#3a6fd0}
#mini{display:block;width:100%;height:56px;background:#0a0c10;border:1px solid #262c38;border-radius:4px;cursor:ew-resize;margin-bottom:4px}
#cv{display:block;width:100%;background:#0a0c10;border:1px solid #262c38;border-radius:4px;cursor:grab}
#cv.grab{cursor:grabbing}
#tip{position:fixed;pointer-events:none;background:#000d;border:1px solid #555;border-radius:4px;padding:6px 8px;font:12px monospace;display:none;z-index:9;white-space:pre;color:#fff}
.legend{display:flex;flex-wrap:wrap;gap:4px 14px;margin:6px 0;font-size:11px;color:#aab}
.legend b{display:inline-block;width:11px;height:11px;margin-right:3px;vertical-align:-1px}
.rtog{cursor:pointer;user-select:none;padding:2px 6px;border-radius:4px;border:1px solid #2a3140}.rtog:hover{background:#ffffff14}.rtog.off{opacity:.32;text-decoration:line-through}
.hint{color:#6b7280;font-size:11px}</style></head><body><div id=wrap>
<h1>Qubit allocation explorer &mdash; __TITLE__</h1>
<div class=bar>
 <button id=mStack class=on onclick="setMode('stack')">Stacked (all roles)</button>
 <button id=mScr onclick="setMode('scratch')">Scratch only</button>
 <span style="width:14px"></span>
 <button onclick="zreset()">Reset</button>
 <button id=pk class=on onclick="togPk()">Peak markers</button>
 <span class=hint id=range></span>
</div>
<canvas id=mini></canvas>
<div class=legend id=leg></div>
<canvas id=cv></canvas>
<div class=hint>Minimap: drag the window to pan, drag empty space to select a span. Main view: scroll/trackpad to zoom (cursor-centered), drag to pan, double-click to reset.</div>
</div><div id=tip></div>
<script>
const DB=__DATA__, R=DB.rows, PH=DB.phases, PEAK=DB.peak, MAXS=DB.maxscr, N=R.length;
let mode='stack', x0=0, x1=N, showPk=true;
const cv=document.getElementById('cv'),x=cv.getContext('2d');
const mini=document.getElementById('mini'),mc=mini.getContext('2d');
const STRIP=22, PAD=26, DH=440, H=DH+PAD+STRIP, W=1500, MH=56;
cv.width=W; cv.height=H; mini.width=W; mini.height=MH;
const ROLES=[['tx',4,'#2b6cd0','tx (dx)'],['ty',5,'#2a9d4f','ty (dy)'],['u',6,'#8e44c0','u (GCD)'],['tr',7,'#e07020','transcript'],['scr',3,'#c0392b','scratch']];
let roleVis={tx:1,ty:1,u:1,tr:1,scr:1};
function phaseOf(p){const s=p.split('/');return s[2]||s[s.length-1]||p;}
function niceStep(top){const raw=top/5;const p=Math.pow(10,Math.floor(Math.log10(raw)));const m=raw/p;const s=m<1.5?1:m<3?2:m<7?5:10;return Math.max(1,s*p);}
function hls(h,l,s){const f=n=>{const k=(n+h*12)%12,a=s*Math.min(l,1-l);return Math.round(255*(l-a*Math.max(-1,Math.min(k-3,Math.min(9-k,1)))));};return[f(0),f(8),f(4)];}
function pcol(p){let h=0;for(const c of p)h=(h*31+c.charCodeAt(0))%360;const[r,g,b]=hls(h/360,0.56,0.55);return`rgb(${r},${g},${b})`;}
function clamp(){x0=Math.max(0,x0);x1=Math.min(N,x1);if(x1-x0<3)x1=x0+3;if(x1>N){x1=N;x0=Math.max(0,N-3);}}
function draw(){
  clamp();
  const n=x1-x0;
  let vmax=1;for(let i=0;i<n;i++){const j=Math.floor(x0)+i;if(j>=N)break;let v;if(mode==='scratch')v=R[j][3];else{v=0;for(const[key,idx]of ROLES)if(roleVis[key])v+=R[j][idx];}if(v>vmax)vmax=v;}
  const top=Math.max(4,Math.ceil(vmax*1.08)), sy=DH/top, yOf=v=>STRIP+PAD+DH-v*sy;
  x.fillStyle='#0a0c10';x.fillRect(0,0,W,H);
  for(let i=0;i<n;i++){const j=Math.floor(x0)+i; if(j>=N)break; const px=i*W/n, pw=Math.ceil(W/n)+1, r=R[j];
    if(mode==='scratch'){x.fillStyle=pcol(phaseOf(PH[r[1]]));x.fillRect(px,yOf(r[3]),pw,DH-(yOf(r[3])-STRIP-PAD));}
    else{let acc=0;for(const[key,idx,col]of ROLES){if(!roleVis[key])continue;const v=r[idx];if(v>0){x.fillStyle=col;x.fillRect(px,yOf(acc+v),pw,v*sy);}acc+=v;}}}
  for(let i=0;i<n;i++){const j=Math.floor(x0)+i;if(j>=N)break;x.fillStyle=pcol(phaseOf(PH[R[j][1]]));x.fillRect(i*W/n,0,Math.ceil(W/n)+1,STRIP-3);}
  x.strokeStyle='#fff2';x.fillStyle='#fff7';x.font='10px monospace';const step=niceStep(top);
  for(let v=step;v<top;v+=step){const y=yOf(v)+.5;x.beginPath();x.moveTo(0,y);x.lineTo(W,y);x.stroke();x.fillText(''+v,3,y-2);}
  x.strokeStyle='#fff8';x.setLineDash([5,4]);x.beginPath();x.moveTo(0,yOf(vmax)+.5);x.lineTo(W,yOf(vmax)+.5);x.stroke();x.setLineDash([]);
  x.fillStyle='#fff';x.fillText((mode==='scratch'?'window max scratch ':'window max active ')+vmax,3,yOf(vmax)-3);
  if(showPk&&mode==='stack'){x.fillStyle='#ff3b30';for(let i=0;i<n;i++){const j=Math.floor(x0)+i;if(j<N&&R[j][2]>=PEAK-1)x.fillRect(i*W/n,STRIP,2,DH+PAD);}}
  const j0=Math.floor(x0),j1=Math.min(N-1,Math.ceil(x1));
  document.getElementById('range').textContent=`snapshots ${j0}–${j1} of ${N}  ·  ops ${R[j0][0].toLocaleString()}–${R[j1][0].toLocaleString()}`;
  if(mode==='stack')document.getElementById('leg').innerHTML=ROLES.map(ro=>`<span class="rtog${roleVis[ro[0]]?'':' off'}" onclick="togRole('${ro[0]}')"><b style="background:${ro[2]}"></b>${ro[3]}</span>`).join('')+'<span style="opacity:.55;margin-left:8px"><b style="background:#ff3b30"></b>peak-touch</span>';
  else{const seen=new Set(),o=[];for(let i=j0;i<=j1;i++){const k=phaseOf(PH[R[i][1]]);if(!seen.has(k)){seen.add(k);o.push(k);}}document.getElementById('leg').innerHTML=o.slice(0,28).map(p=>`<span><b style="background:${pcol(p)}"></b>${p}</span>`).join('');}
  drawMini();
}
function drawMini(){
  mc.fillStyle='#0a0c10';mc.fillRect(0,0,W,MH);
  mc.fillStyle='#34507e';for(let j=0;j<N;j++){const h=R[j][2]/PEAK*(MH-2);mc.fillRect(j*W/N,MH-h,Math.ceil(W/N)+1,h);}
  const a=x0/N*W,b=x1/N*W;mc.fillStyle='#ffffff26';mc.fillRect(a,0,b-a,MH);
  mc.strokeStyle='#5b8ff0';mc.lineWidth=2;mc.strokeRect(a+1,1,Math.max(2,b-a-2),MH-2);
}
function setMode(m){mode=m;document.getElementById('mStack').className=m==='stack'?'on':'';document.getElementById('mScr').className=m==='scratch'?'on':'';draw();}
function zreset(){x0=0;x1=N;draw();}
function setRange(a,b){x0=a;x1=b;clamp();draw();}
function togPk(){showPk=!showPk;document.getElementById('pk').className=showPk?'on':'';draw();}
function togRole(k){roleVis[k]=!roleVis[k];draw();}
// wheel/trackpad zoom on main, cursor-centered
cv.addEventListener('wheel',e=>{e.preventDefault();const r=cv.getBoundingClientRect();const frac=Math.max(0,Math.min(1,(e.clientX-r.left)/r.width));const cx=x0+frac*(x1-x0);let f=Math.exp(e.deltaY*0.0016);f=Math.max(0.5,Math.min(2,f));x0=cx-(cx-x0)*f;x1=cx+(x1-cx)*f;draw();},{passive:false});
// drag-pan on main
let pan=null;
cv.addEventListener('mousedown',e=>{const r=cv.getBoundingClientRect();pan={f:(e.clientX-r.left)/r.width,x0,x1};cv.classList.add('grab');});
cv.addEventListener('dblclick',zreset);
// minimap: drag window to pan, drag empty to brush
let md=null;
mini.addEventListener('mousedown',e=>{const r=mini.getBoundingClientRect();const frac=(e.clientX-r.left)/r.width,px=frac*N;if(px>=x0&&px<=x1)md={m:'pan',f:frac,x0,x1};else md={m:'brush',s:px};});
window.addEventListener('mousemove',e=>{
  if(pan){const r=cv.getBoundingClientRect();const frac=(e.clientX-r.left)/r.width;const d=(pan.f-frac)*(pan.x1-pan.x0);x0=pan.x0+d;x1=pan.x1+d;clamp();draw();return;}
  if(md){const r=mini.getBoundingClientRect();const frac=(e.clientX-r.left)/r.width,px=frac*N;if(md.m==='pan'){const d=(frac-md.f)*N;x0=md.x0+d;x1=md.x1+d;}else{x0=Math.min(md.s,px);x1=Math.max(md.s,px);}clamp();draw();return;}
});
window.addEventListener('mouseup',()=>{pan=null;md=null;cv.classList.remove('grab');});
// hover tooltip (only when not dragging)
const tip=document.getElementById('tip');
cv.addEventListener('mousemove',e=>{if(pan)return;const r=cv.getBoundingClientRect();const i=Math.floor((e.clientX-r.left)/r.width*(x1-x0));const j=Math.floor(x0)+i;if(j<0||j>=N){tip.style.display='none';return;}const d=R[j];
  tip.textContent=`op ~${d[0].toLocaleString()}\nactive ${d[2]}\nscratch ${d[3]} | tx ${d[4]} ty ${d[5]} u ${d[6]} tr ${d[7]}\nphase: ${phaseOf(PH[d[1]])}`;
  tip.style.display='block';tip.style.left=(e.clientX+14)+'px';tip.style.top=(e.clientY+14)+'px';});
cv.addEventListener('mouseleave',()=>tip.style.display='none');
draw();
</script></body></html>"""
open(OUT, "w").write(PAGE.replace("__DATA__", D).replace("__TITLE__", TITLE))
print(f"wrote {OUT}  (N={len(rows)}, peak={peak}, maxscr={maxscr})")
