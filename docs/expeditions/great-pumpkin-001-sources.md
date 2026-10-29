# Sources and comparison with published work

Sources accessed 2026-10-28. Dates below are publication/revision dates, not
search-engine crawl dates. This is a focused primary-source review, not an
exhaustive systematic review or a novelty claim.

| Source | What was checked / implication |
|---|---|
| [Luo, Wang, Xiao, *Two-loop Renormalization Group Equations in General Gauge Field Theories*, PRD 67 (2003) 065019; hep-ph/0211440v3](https://arxiv.org/pdf/hep-ph/0211440) | Eq. (30), gauge beta function, and product-group substitution rules. Real-scalar conventions convert to the code's complex-scalar weights. Yukawa terms are explicitly present in the literature and absent from the scanner. |
| [Bertolini, Di Luzio, Malinský, *Intermediate mass scales in the non-supersymmetric SO(10) grand unification: a reappraisal*, PRD 80 (2009) 015013; 0903.4049v2, 16 June 2009](https://arxiv.org/pdf/0903.4049) | Table VI, row Ia and other identical 3221 entries: exact one/two-loop gauge matrix after permuting `(L,R,BL,3)` to `(3,L,R,BL)`. Table VIII: SM coefficients. Eq. (26): finite one-loop matching constants. Sections III D/E discuss small Yukawa corrections in their examples and sensitivity to thresholds. The expedition's approximate `log10 MG≈16.2`, `1/αG≈45.4` are familiar, not novel. Different input/matching and intermediate-stage assumptions prevent a precision replication claim for their scale plots. |
| [Bertolini, Di Luzio, Malinský, *Seesaw Scale in the Minimal Renormalizable SO(10) Grand Unification*, PRD 85 (2012) 095014; 1202.0807](https://arxiv.org/abs/1202.0807) | A 45+126 scalar-potential analysis finds higher B−L scales with non-minimal-survival spectra. This prevents interpreting our factor-three, fixed-light-content obstruction as a general SO(10) no-go. We did not reproduce their scalar-vacuum analysis. |
| [Jarkovská, Malinský, Mede, Susič, *Quantum nature of the minimal potentially realistic SO(10) Higgs model*, PRD 105 (2022) 095003; 2109.06784v2, 27 August 2022](https://arxiv.org/abs/2109.06784) | Full scalar corrections and perturbativity constrain the 45+126 model much more strongly than arbitrary mass boxes; the 421 route is preferred in their study. Our independently varied scalar masses are a relaxed domain, not a set of certified vacua. |
| [Jarkovská, Malinský, Susič, *The trouble with the minimal renormalizable SO(10) GUT*, PRD 108 (2023) 055003; 2304.14227](https://arxiv.org/pdf/2304.14227) | Reports difficulty accommodating the required light SM Higgs doublet in perturbative 45+126+10C. This is a model-consistency challenge even when gauge trajectories meet. Our calculation does not resolve or refute it. |
| [Preda, Senjanović, Zantedeschi, *SO(10) theory on the plateau: the importance of being renormalizable*, 2502.21180v1, 28 February 2025](https://arxiv.org/pdf/2502.21180) | Studies additional light scalar spectra and low LR scales; discusses the preceding Higgs issue and calls for independent confirmation/remedies. The broader spectra and assumptions differ from our minimal-survival box. No blanket claim that the model is either established or definitively excluded is warranted here. |
| [PDG, *Electroweak Model and Constraints on New Physics*, revised November 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf) | Section 10.3.1: MS electromagnetic inverse coupling 127.930±0.008; Table 10.5: MS weak angle 0.23122±0.00006; Table 10.6: αs=0.1177±0.0009 and MZ=91.1880 GeV. These are a scheme-aware input sensitivity comparison, not a combined independent Gaussian likelihood. The electromagnetic value carries the review's active-flavour conventions. A complete electroweak matching treatment is not supplied. |
| [Super-Kamiokande, *Search for proton decay via p→e+π0 and p→μ+π0 with an enlarged fiducial volume*, PRD 102 (2020) 112011; 2010.16098v2](https://arxiv.org/abs/2010.16098) | Channel-specific 90% CL partial-lifetime limits 2.4×10^34 years and 1.6×10^34 years. Use the former for e+π0; the repository's historical 1.6×10^34 e+π0 value is stale. |
| [Tashiro for Super-Kamiokande, Moriond EW 2025, slide 14](https://indico.in2p3.fr/event/35965/contributions/152404/attachments/91535/139517/9_TTashiro-v1.pdf) | Collaboration presentation still quotes 2.4×10^34 years for e+π0. The focused search found no stronger channel-specific primary limit to use; we do not claim exhaustive coverage of all 2026 proceedings. |
| [Lean 4 official reference: linear arithmetic solver](https://lean-lang.org/doc/reference/latest/The--grind--tactic/Linear-Arithmetic-Solver/) and [4.24.0 release](https://github.com/leanprover/lean4/releases/tag/v4.24.0) | Standard ordered-field algebra is sufficient for this experiment. The actual proof pins 4.24.0, so it does not depend on future changes to the `latest` reference. |

## What is and is not new here

Gauge coupling differences, linear one-loop scale solving, universal complete-
multiplet shifts, threshold sensitivity and the general success/failure patterns
of non-SUSY SO(10) have long been studied. The expedition does not name these
as discoveries of new physics.

The new repository-level evidence is: exact reproduction and diagnosis of the
historical projected-residual degeneracy; explicit numerical comparison to an
independent inverse-coupling solver; a documented 16-scalar threshold domain;
its `5.4 log10 f` support bound versus the correlated `21/145 log10 f`
response; and a small machine-checked algebraic explanation. Whether that
particular packaging or spectrum-bound presentation is publication-novel was
not established and is not claimed.

The main scientific comparison is a limitation: published scalar-potential
studies address a requirement our gauge-sector calculation cannot satisfy by
itself. Matching their physically allowed mass correlations is the next
investigation, not another unconstrained residual scan.
