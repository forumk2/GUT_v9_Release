"""Static scientific figure from saved evidence; no live fitting."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

HERE=Path(__file__).resolve().parent
d=json.loads((HERE/"results/experiments.json").read_text())
t=d["thresholds"]
box=t["factor3_exact_one_loop_box"]
response=np.array(box["response"])
rho=np.log10(3)
# Exact affine image of the 16-dimensional box: build the 2D zonotope by
# Minkowski addition, retaining convex hull vertices at each step.
vertices=np.array([[0.,0.]])
for col in response.T:
    vertices=np.vstack([vertices+rho*col,vertices-rho*col])
    if len(vertices)>=3 and np.linalg.matrix_rank(vertices-vertices[0])==2:
        vertices=vertices[ConvexHull(vertices).vertices]
vertices+=np.array(box["center"])
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout="constrained")
ax=axs[0]
ax.fill(vertices[:,0],vertices[:,1],alpha=.22,color="#3b82a0",label="Exact one-loop mass-box image")
pts=np.array([[x["result"]["log10_MI"],x["result"]["log10_MGUT"]] for x in t["samples"]])
ax.scatter(*pts.T,s=14,color="#b95127",alpha=.75,label="64 two-loop samples")
for label,color,key in [("Degenerate (one loop)","#174967","one_loop"),("Degenerate (two loop)","#8b321a","two_loop")]:
    row=d[key]; ax.scatter(row["log10_MI"],row["log10_MGUT"],marker="*",s=120,color=color,label=label)
corner=t["local_two_loop_max"]["result"]
ax.scatter(corner["log10_MI"],corner["log10_MGUT"],marker="x",s=65,c="#8b321a",label="Two-loop local maximum of MI")
ax.set(xlabel=r"$\log_{10}(M_I/\mathrm{GeV})$",ylabel=r"$\log_{10}(M_G/\mathrm{GeV})$",
       title="Independent heavy masses within a factor of 3")
ax.legend(fontsize=8,loc="upper right")
ax=axs[1]
factors=np.linspace(1,10,250); center=d["one_loop"]["log10_MI"]
half=5.4*np.log10(factors)
common=abs(sum(response[0]))*np.log10(factors)
ax.fill_between(factors,center-half,center+half,color="#3b82a0",alpha=.22,label="Independent heavy-mass box")
ax.fill_between(factors,center-common,center+common,color="#174967",alpha=.7,label="One common heavy scalar mass")
ax.axhline(13,color="#b95127",linestyle="--",label="Illustrative target MI = 10¹³ GeV")
ax.axvline(3,color="gray",linestyle=":")
ax.set(xlabel="Allowed mass factor f: MG/f ≤ MH ≤ f MG",ylabel=r"One-loop $\log_{10}(M_I/\mathrm{GeV})$",
       title="Mass correlations control predictive power",ylim=(3,15))
ax.legend(fontsize=8,loc="lower left")
fig.suptitle("Expedition 001 · Conditional gauge-sector results",fontsize=13)
fig.savefig(HERE/"results/threshold_geometry.png",dpi=180)
fig.savefig(HERE/"results/threshold_geometry.svg")
svg=HERE/"results/threshold_geometry.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n",
               encoding="utf-8",newline="\n")
print("wrote threshold_geometry.png and .svg")
