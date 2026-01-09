## Note on Numerical Differences vs. the Paper
## What Changed (Important Implementation Note)

According to the collision-model formulation used in the paper, when resolving the momentum transfer on a candidate impact face, the following three vectors must be **coplanar**:

- **Relative velocity vector**: `v_r`
- **Impact-face unit normal**: `n`
- **Impact-face unit tangential direction**: `t`

Mathematically, this requirement can be written as:

v_r ∈ span{n, t}
⇔
n · (t × v_r) = 0


In the original implementation (used for the published results), the tangential direction `t` was constructed using a **fixed global reference axis**, which does not generally guarantee that `v_r`, `n`, and `t` lie in the same plane. This issue is particularly noticeable in **BIP (eccentric impact)** cases.

The current code fixes this by constructing `t` from the **incoming relative velocity direction**, i.e. the projection of `-v_r` onto the impact-face tangent plane. This enforces coplanarity **by construction**.

### Impact on the Results

This issue only affected the **quantitative values** of the computed deflection.  
The **qualitative conclusions of the paper remain valid**, in particular:

**Eccentric (off-center) impacts can provide additional deflection benefits compared to centric impacts.**

After correcting this geometric inconsistency, the updated results are **more geometrically consistent** and may yield **slightly improved** deflection performance.

### Example: (99942) Apophis

As a concrete example, the asteroid **(99942) Apophis** highlights the numerical impact of the corrected implementation.

Using the original (paper-version) code:

- **COG strategy deflection distance**: 1,306,060.18  
- **BIP strategy deflection distance**: 1,769,910.92  
- **Relative gain (BIP vs. COG)**: **35.52%**

After fixing the coplanarity issue between `v_r`, `n`, and `t`:

- **COG strategy deflection distance**: 1,359,482.65  
- **BIP strategy deflection distance**: 1,957,621.50  
- **Relative gain (BIP vs. COG)**: **44%**

These results show that enforcing the correct geometric relationship does not change the qualitative conclusion of the paper, but leads to **quantitatively larger and more physically consistent deflection estimates**, particularly for the eccentric (BIP) impact strategy.


