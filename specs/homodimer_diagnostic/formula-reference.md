# Verified Formula Reference — 7 Complex Confidence Values

**Ground truth:** `specs/homodimer_diagnostic/references/ipsae_v4.py`
(DunbrackLab/IPSAE v4, 3 Jan 2026, md5 `a48df7adc64afa36b4b475429b0a8aff`, 1010 lines).
All line citations below are `ipsae_v4.py:<line>` and were read directly from that file.

**Compared against:**
- `notebooks/homodimer_diagnostic.ipynb` cell 18 (`# Core scoring functions`) — "notebook" below.
- `CLAUDE.md` §"Skill: Homodimer Confidence Scoring" — "CLAUDE.md" below.

**Scope note.** `ipsae.py` emits three rows per chain pair: two `asym` rows (A→B and B→A) and
one `max` row. The single number a user quotes is the `max` row (`ipsae_v4.py:983-1002`).
Every "authoritative" formula below is the `max`-row value unless stated otherwise.

**Conventions.** `pae_AB = pae_matrix[:nA, nA:]` (rows = chain A, cols = chain B);
`pae_BA = pae_matrix[nA:, :nA]`. PAE is asymmetric, so `pae_AB != pae_BA.T`.
Default `pae_cutoff = 10.0` (CLI arg, `ipsae_v4.py:58`), contact cutoff `8.0 Å`
(`pDockQ_cutoff`, `ipsae_v4.py:644`), LIS cutoff hardcoded `12` (`ipsae_v4.py:712`).

---

## 0. The two `d0` helpers — NOT equivalent

```python
# ipsae_v4.py:116-124
def calc_d0(L, pair_type):            # scalar
    L = float(L)
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    d0 = 1.24 * (L - 15) ** (1/3) - 1.8 if L > 27 else 1.0
    return max(min_value, d0)

# ipsae_v4.py:126-136
def calc_d0_array(L, pair_type):      # vectorised
    L = np.maximum(26, np.array(L, dtype=float))     # clamp, ipsae_v4.py:130
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    return np.maximum(min_value, 1.24 * (L - 15) ** (1/3) - 1.8)
```

`ptm_func(x, d0) = 1.0 / (1 + (x/d0)**2)` (`ipsae_v4.py:111-112`); vectorised via
`np.vectorize` as `ptm_func_vec` (`ipsae_v4.py:113`).

**Which helper is used where:**

| d0 | helper | citation |
|----|--------|----------|
| `d0chn` (ipTM_d0chn, ipSAE_d0chn) | `calc_d0` scalar | `ipsae_v4.py:732` |
| `d0dom` (ipSAE_d0dom) | `calc_d0` scalar | `ipsae_v4.py:778` |
| `d0res` (ipSAE_d0res) | `calc_d0_array` vectorised | `ipsae_v4.py:787` |
| pDockQ2 | neither — fixed `d0 = 10.0` | `ipsae_v4.py:687` |

**Divergence.** Enumerated for all integer `L` in 0..2000, the two helpers agree **exactly**
everywhere except `L == 27`:

| L | `calc_d0` | `calc_d0_array` |
|---|-----------|-----------------|
| ≤ 26 | 1.0 | 1.0 (clamp → 26 → cubic = 0.9578 → floored to 1.0) |
| **27** | **1.0** (the `L > 27` test at `ipsae_v4.py:120` is false) | **1.038891** |
| ≥ 28 | cubic | cubic (clamp inactive) |

So the clamp at `ipsae_v4.py:130` is not what creates the difference — the `max(1.0, …)`
floor already handles `L < 26`. The difference is created entirely by the `else: d0 = 1.0`
branch at `ipsae_v4.py:122-123`, which `calc_d0_array` has no counterpart for. `L == 27` is
the only point where the cubic exceeds 1.0 while the scalar branch still returns 1.0.

- **notebook:** defines one `d0_func` with the `L <= 27` branch (scalar semantics) and uses it
  for `d0chn`, `d0dom`, **and** `d0res`. Correct for the first two, wrong for `d0res` when a
  residue has exactly 27 valid inter-chain pairs.
- **CLAUDE.md:** defines only the scalar `d0_func`; same defect, and does not mention that a
  second helper exists.
- **Verdict:** ipsae.py is authoritative. Both derivatives need a separate `d0_array` for
  `d0res` only. The `pair_type='nucleic_acid'` → `min_value = 2.0` branch
  (`ipsae_v4.py:119, 133`) is inert for protein–protein pairs; neither derivative models it,
  which is acceptable for this project's homodimer scope.

---

## 1. ipTM_d0chn

Per-residue mean of `ptm(PAE, d0chn)` over **all** residues of the partner chain — **no PAE
cutoff**. `d0chn = calc_d0(nA + nB)`.

```python
# ipsae_v4.py:731-732, 734, 736, 746, 827-829, 851
n0chn = nA + nB                                    # :731
d0chn = calc_d0(n0chn, 'protein')                  # :732  scalar helper
ptm_d0chn = ptm_func(pae_matrix, d0chn)            # :734
valid_pairs_iptm = (chains == chain2)              # :736  ALL partner residues, no cutoff
iptm_byres[i] = ptm_d0chn[i, valid_pairs_iptm].mean()   # :746
iptm_asym[A][B] = iptm_byres.max()                 # :827-829 argmax over rows of chain1
iptm_max = max(iptm_asym[A][B], iptm_asym[B][A])   # :851
```

| source | behaviour |
|--------|-----------|
| ipsae.py | as above; `d0chn` from the **scalar** helper (`:732`) |
| notebook | identical: `d0_func(nA+nB)`, `ptm_func(pae_ab, d0).mean(axis=1).max()`, then `max` over directions |
| CLAUDE.md | identical formula |

**Verdict: all three agree. Naming only** — this is `ipTM_d0chn`, a reimplementation from the
PAE matrix. It is **not** AlphaFold's own `ipTM` (`iptm_af`, `ipsae_v4.py:944-946`), which is
read from the model's summary file and which AFDB does not expose on these endpoints. The
notebook labels it plain `ipTM` (cell 19); that label is misleading.

---

## 2. ipSAE_d0chn

Same as ipTM_d0chn but restricted to pairs with `PAE < pae_cutoff`. Same `d0chn`.

```python
# ipsae_v4.py:737, 747, 832-834, 859
valid_pairs_matrix = np.outer(chains == chain1, chains == chain2) & (pae_matrix < pae_cutoff)  # :737
ipsae_d0chn_byres[i] = ptm_d0chn[i, valid_pairs_matrix[i]].mean() if any else 0.0              # :747
ipsae_d0chn_asym = byres.max()                                                                  # :832-834
ipsae_d0chn_max  = max(asym[A][B], asym[B][A])                                                  # :859
```

Note the cutoff test is strict `<` (`ipsae_v4.py:737`), and `d0chn` uses `calc_d0`
(`ipsae_v4.py:732`), not the array helper.

| source | behaviour |
|--------|-----------|
| ipsae.py | as above |
| notebook | identical (`r_d0chn[i] = ptm_func(vp, d0_chn).mean()`, strict `<`, max over directions) |
| CLAUDE.md | identical formula (`d0_mode == 'chn'` branch) |

**Verdict: all three agree.**

---

## 3. ipSAE_d0res

Per row `i`, `d0` is derived from that row's own count of valid pairs, using the **array**
helper.

```python
# ipsae_v4.py:783, 786-787, 795, 798-800, 842-844, 885
valid_pairs_matrix = np.outer(chains == chain1, chains == chain2) & (pae_matrix < pae_cutoff)  # :783
n0res_byres = valid_pairs_matrix.sum(axis=1)                    # :786
d0res_byres = calc_d0_array(n0res_byres, pair_type)             # :787  ARRAY helper
ptm_row     = ptm_func(pae_matrix[i], d0res_byres[i])           # :798-799
ipsae_d0res_byres[i] = ptm_row[valid_pairs_matrix[i]].mean() if any else 0.0   # :800
ipsae_d0res_asym = byres.max()                                  # :842-844
ipsae_d0res_max  = max(asym[A][B], asym[B][A])                  # :885
```

The reported `n0res` / `d0res` are those of the argmax residue (`ipsae_v4.py:846-847`), and
for the `max` row those of the winning direction (`ipsae_v4.py:888-893, 898-901`).

| source | behaviour |
|--------|-----------|
| ipsae.py | `calc_d0_array(n0res)` per residue |
| notebook | `d0_func(n0)` — the **scalar** helper |
| CLAUDE.md | `d0_func(n0)` — the **scalar** helper |

**Verdict: ipsae.py is authoritative.** Notebook and CLAUDE.md are correct for every residue
except one with exactly 27 valid inter-chain pairs, where they use `d0 = 1.0` instead of
`1.038891`. Because a residue with only 27 sub-cutoff partners is rarely the argmax, the
reported score is usually unaffected; but it can differ, and the project tolerance is ±0.001.
Note the magnitude 0.0389 is the **d0** difference, not the score difference — the induced
score difference is smaller and depends on the PAE distribution in that row.

---

## 4. ipSAE_d0dom

`d0dom` is derived from the number of residues (both chains) participating in at least one
sub-cutoff inter-chain pair — computed **per direction** from that direction's own PAE block.

```python
# ipsae_v4.py:737, 749-756, 775-778, 781, 796, 837-839, 867
# sets accumulated during the chain1->chain2 pass, from valid_pairs_matrix (:737):
#   unique_residues_chain1[A][B] = rows of the A->B block with >= 1 valid pair   # :751-753
#   unique_residues_chain2[A][B] = cols of the A->B block with >= 1 valid pair   # :754-756
n0dom = len(unique_residues_chain1[A][B]) + len(unique_residues_chain2[A][B])    # :775-777
d0dom = calc_d0(n0dom, pair_type)                                               # :778  scalar helper
ptm_d0dom = ptm_func(pae_matrix, d0dom)                                         # :781
ipsae_d0dom_byres[i] = ptm_d0dom[i, valid_pairs_matrix[i]].mean() if any else 0.0  # :796
ipsae_d0dom_asym = byres.max()                                                  # :837-839
ipsae_d0dom_max  = max(asym[A][B], asym[B][A])                                  # :867
```

Equivalent block form for direction X→Y with block `P`:

```python
n0dom_XY = int((P < pae_cutoff).any(axis=1).sum()) + int((P < pae_cutoff).any(axis=0).sum())
```

Because PAE is asymmetric, `n0dom(A→B) != n0dom(B→A)` in general, so `d0dom` differs between
the two `asym` rows. The `max` row reports the `n0dom`/`d0dom` of the winning direction
(`ipsae_v4.py:870-875, 880-883`).

| source | behaviour |
|--------|-----------|
| ipsae.py | per direction: rows-with-any + cols-with-any of **that direction's own block** |
| notebook | computes `n_dom` **once** from `pae_ab` (`any(axis=1).sum() + any(axis=0).sum()`) and passes the same `d0_dom` into both `_per_res(pae_ab, nA)` and `_per_res(pae_ba, nB)` |
| CLAUDE.md | `n_dom_A = (pae_AB < c).any(axis=1).sum()`, `n_dom_B = (pae_BA < c).any(axis=1).sum()`, `n_dom = max(n_dom_A + n_dom_B, 1)`, one value for both directions |

**Verdict: ipsae.py is authoritative; both derivatives are wrong, differently.**
- Notebook: correct value for the A→B direction, reused for B→A. The B→A per-residue array is
  therefore computed with the wrong `d0`.
- CLAUDE.md: mixes *rows of `pae_AB`* with *rows of `pae_BA`*. The second term should be
  *columns of `pae_AB`* for the A→B direction. Both terms happen to count "residues of B with
  some good inter-chain PAE", but measured in different blocks, so the counts differ. Also
  single-valued across directions, same defect as the notebook. The `max(…, 1)` clamp has no
  counterpart in ipsae.py; it is harmless (`calc_d0(0)` already returns 1.0).

---

## 5. pDockQ

`x` is built from the **number of contact pairs**, not the number of interface residues.

```python
# ipsae_v4.py:643-670
pDockQ_cutoff = 8.0                                             # :644
npairs = 0
for i in rows_of_chain1:
    valid_pairs = (chains == chain2) & (distances[i] <= 8.0)    # :652  note '<=' 
    npairs += valid_pairs.sum()                                 # :653  PAIR count
    if valid_pairs.any():
        unique.add(i); unique.update(np.where(valid_pairs)[0])  # :654-659
if npairs > 0:
    mean_plddt = cb_plddt[list(unique)].mean()                  # :663  union over BOTH chains
    x = mean_plddt * math.log10(npairs)                         # :664
    pDockQ = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018   # :665
else:
    pDockQ = 0.0                                                # :669
```

`nres = len(unique)` is computed at `ipsae_v4.py:662` but is **only printed**
(`ipsae_v4.py:963-964`); it never enters the score.

`distances` are CB–CB, with CA substituted for GLY (`ipsae_v4.py:385` selects
`CB`, or `C3*` for nucleic acids, or `CA` when `residue_name == "GLY"`;
`ipsae_v4.py:402, 425`). `cb_plddt` is the pLDDT of that CB atom (`ipsae_v4.py:534`).

pDockQ is **direction-independent**: `npairs` and the residue set are both symmetric under
swapping the chains, which is why the `max` row prints `pDockQ[chain1][chain2]` unchanged
(`ipsae_v4.py:989`) rather than a max.

| source | behaviour |
|--------|-----------|
| ipsae.py | `x = mean_plddt * log10(npairs)`, `npairs` = contact **pairs** |
| notebook | `n_c = int(cmask.sum())` = contact pairs; `x = mean_plddt * np.log10(n_c)`; mean pLDDT over `plddt_a[if_A]` ∪ `plddt_b[if_B]`; same constants — **matches** |
| CLAUDE.md | `n = int(if_A.sum()) + int(if_B.sum())` = interface **residue count**; `x = mean_plddt * np.log10(max(n, 1))` — **wrong**, this is ipsae.py's unused `nres` |

**Verdict: notebook correct (except the zero-contact return, §8); CLAUDE.md wrong.** The
draft's claim that pDockQ uses contact pairs is **confirmed** at `ipsae_v4.py:653, 664`.
For a typical interface `npairs` is several times `nres`, so CLAUDE.md's `x` is materially
too small.

---

## 6. pDockQ2

```python
# ipsae_v4.py:672-700
npairs = 0; total = 0.0
for i in rows_of_chain1:
    valid_pairs = (chains == chain2) & (distances[i] <= 8.0)    # :683
    if valid_pairs.any():
        npairs += valid_pairs.sum()                             # :685
        total  += ptm_func(pae_matrix[i][valid_pairs], 10.0).sum()   # :686-688  fixed d0 = 10.0
if npairs > 0:
    mean_plddt = cb_plddt[list(pDockQ_unique_residues[A][B])].mean()  # :692  same set as pDockQ
    mean_ptm   = total / npairs                                 # :693
    x = mean_plddt * mean_ptm                                   # :694
    pDockQ2 = 1.31 / (1 + math.exp(-0.075 * (x - 84.733))) + 0.005    # :695
else:
    pDockQ2 = 0.0                                               # :700
```

Constants **confirmed** at `ipsae_v4.py:695`: numerator `1.31`, slope `-0.075`, midpoint
`84.733`, offset `+0.005`.

pDockQ2 **is** directional — `mean_ptm` reads only the X→Y PAE block (`ipsae_v4.py:686`),
while `mean_plddt` is symmetric. The `max` row therefore reports
`max(pDockQ2[A][B], pDockQ2[B][A])` (`ipsae_v4.py:977, 990`).

`pDockQ_unique_residues` is populated in the pDockQ loop (`ipsae_v4.py:655-659`) and reused
here; the pDockQ block must run first.

| source | behaviour |
|--------|-----------|
| ipsae.py | `mean_ptm` over contact pairs at `d0 = 10`, directional; `mean_plddt` over unique interface residues; `1.31 / (1 + exp(-0.075(x - 84.733))) + 0.005`; reported value = max over directions |
| notebook | correct constants and correct `mean_plddt`; but `compute_pdockq2(dist_matrix, pae_AB, …)` (cell 19) evaluates the **A→B direction only** and reports it as the final score |
| CLAUDE.md | `0.715 / (1 + exp(-12.3 * (x - 0.605))) + 0.005` — **wrong constants**; and `plddt_contacts` is built with `np.repeat` / `np.tile` over the contact mask, i.e. pLDDT weighted by each residue's contact multiplicity rather than a mean over unique interface residues — **also wrong** |

**Verdict: ipsae.py authoritative. Notebook has the right formula but reports the asymmetric
A→B value instead of `max(A→B, B→A)`. CLAUDE.md is wrong on both the sigmoid constants and
the pLDDT aggregation.**

---

## 7. LIS

```python
# ipsae_v4.py:702-720
mask = (chains[:, None] == chain1) & (chains[None, :] == chain2)   # :708  directional block
selected_pae = pae_matrix[mask]                                     # :709
valid_pae = selected_pae[selected_pae < 12]                         # :712  hardcoded 12, strict '<'
LIS[A][B] = ((12 - valid_pae) / 12).mean() if valid_pae.size else 0.0   # :713-718

# reported score:
LIS_Score = (LIS[A][B] + LIS[B][A]) / 2.0                           # :982  MEAN, not max
```

The cutoff `12` is hardcoded at `ipsae_v4.py:712, 714` and is independent of the `pae_cutoff`
CLI argument. The `asym` rows print the directional values (`ipsae_v4.py:956`); the `max` row
prints the mean (`ipsae_v4.py:991`).

| source | behaviour |
|--------|-----------|
| ipsae.py | per-direction `mean((12 - PAE)/12)` over `PAE < 12`; reported = **mean** of the two directions |
| notebook | `compute_lis` returns `(lis_ab + lis_ba) / 2.0` — **matches** |
| CLAUDE.md | `score = max(_lis_half(pae_AB), _lis_half(pae_BA))` — **wrong**, uses max |

**Verdict: notebook correct; CLAUDE.md wrong.** The draft's claim that LIS combines by mean
is **confirmed** at `ipsae_v4.py:982`. `max ≥ mean` always, so CLAUDE.md systematically
overestimates LIS.

---

## 8. Zero inter-chain contacts / zero valid pairs

| condition | ipsae.py | notebook | CLAUDE.md |
|-----------|----------|----------|-----------|
| no CB–CB contact ≤ 8 Å | `pDockQ = 0.0` (`ipsae_v4.py:666-669`) | `0.018` | `0.018` |
| no CB–CB contact ≤ 8 Å | `pDockQ2 = 0.0` (`ipsae_v4.py:696-700`) | `0.005` | `0.005` |
| no inter-chain `PAE < 12` | `LIS = 0.0` (`ipsae_v4.py:717-720`) | `0.0` | `0.0` |
| residue with no `PAE < cutoff` pair | that residue's ipSAE entry stays `0.0` (`ipsae_v4.py:747, 796, 800`) | `0.0` | `0.0` |
| ipTM_d0chn | never degenerate — `valid_pairs_iptm = (chains == chain2)` (`ipsae_v4.py:736`) is non-empty for any real chain pair, so the per-residue mean is always defined |

**Verdict: ipsae.py returns `0.0`, not the sigmoid minimum.** `0.018` and `0.005` are the
`x → -∞` limits of the published Bryant and Zhu sigmoids and are what those authors' own code
returns; `ipsae.py` deliberately short-circuits to `0.0` instead. The notebook and CLAUDE.md
follow the published-sigmoid convention. This is a **defensible disagreement**, not a
transcription error — but it must be chosen once and documented, because the two conventions
give different answers for a non-interacting pair. If ±0.001 agreement with `ipsae.py` is the
acceptance criterion, `0.0` is the required behaviour.

---

## Summary table

| Value | ipsae.py citation | notebook | CLAUDE.md |
|-------|-------------------|----------|-----------|
| ipTM_d0chn | `:732, :736, :746, :851` | correct (mislabelled `ipTM`) | correct |
| ipSAE_d0chn | `:732, :737, :747, :859` | correct | correct |
| ipSAE_d0res | `:786-787, :798-800, :885` | wrong `d0` helper at `n0res == 27` | wrong `d0` helper at `n0res == 27` |
| ipSAE_d0dom | `:749-756, :775-778, :796, :867` | `n_dom` not recomputed for B→A | rows of `pae_AB` + rows of `pae_BA`, and not per direction |
| pDockQ | `:652-653, :663-665` | correct | uses residue count instead of pair count |
| pDockQ2 | `:683-695, :977` | correct formula, reports A→B not `max` | wrong constants + wrong pLDDT aggregation |
| LIS | `:708-718, :982` | correct | uses `max` instead of `mean` |
| `calc_d0` vs `calc_d0_array` | `:116-124` vs `:126-136` | only the scalar form exists | only the scalar form exists |
| zero contacts | `:669, :700, :718` | sigmoid minima | sigmoid minima |

---

## Corrections to the draft audit

Verified row by row against `ipsae_v4.py` source. The draft in
`specs/homodimer_diagnostic/rework-plan.md` §R001 is **substantively correct on every
formula claim**, including all five claims flagged for particular attention:

1. pDockQ uses the count of contact **pairs** — **confirmed** (`ipsae_v4.py:653, 664`).
   `nres` is computed (`:662`) but never used in the score.
2. pDockQ2 constants `1.31 / -0.075 / 84.733` — **confirmed** (`ipsae_v4.py:695`).
3. LIS combines as the **mean** of the two directions — **confirmed** (`ipsae_v4.py:982`).
4. `n0dom` is per direction, rows-with-any + cols-with-any of that same block — **confirmed**
   (`ipsae_v4.py:751-756, 775-777`).
5. `calc_d0_array` clamps `L` to ≥ 26 while `calc_d0` has an `L <= 27` branch, differing only
   at `L == 27` (1.038891 vs 1.0) — **confirmed** (`ipsae_v4.py:120-123` vs `:130`).
6. `ipsae.py` returns `0.0`, not the sigmoid minimum, when there are no contacts —
   **confirmed** (`ipsae_v4.py:669, 700`).

The following defects in the draft were found and are corrected above:

**D1 — One substantive omission: pDockQ2 direction handling.** The draft's pDockQ2 row says
"notebook correct". It is not fully correct. `ipsae.py` computes pDockQ2 per direction
(`mean_ptm` reads only the X→Y PAE block, `ipsae_v4.py:686`) and reports
`max(pDockQ2[A][B], pDockQ2[B][A])` (`ipsae_v4.py:977, 990`). The notebook calls
`compute_pdockq2(dist_matrix, pae_AB, …)` (cell 19) and reports the A→B value only. This is a
real numerical difference of the same kind as the `d0dom` bug the draft *did* catch, and it
has no corresponding task in the plan (R003 covers `d0dom`, nothing covers pDockQ2). The
draft's own phrase "using the directional block" describes ipsae.py accurately but does not
notice the notebook mismatch that follows from it.

**D2 — Line citations are off by one for ipTM_d0chn and ipSAE_d0chn.** The draft cites
`L747` for ipTM_d0chn and `L748` for ipSAE_d0chn. The actual lines are **`:746`**
(`iptm_d0chn_byres`) and **`:747`** (`ipsae_d0chn_byres`).

**D3 — Incomplete citation for ipSAE_d0dom.** The draft cites `L775-778`, which is only where
`n0dom` is summed and `d0dom` computed. The sets it sums are built much earlier, in the ipTM
loop at **`:749-756`**, from `valid_pairs_matrix` (`:737`). Anyone reimplementing from
`L775-778` alone cannot see what is being counted.

**D4 — Minor: `calc_d0_array` description is incomplete.** The draft says it "clamps
`L = max(26, L)` first, then always applies the cubic form". It also applies a
`max(min_value, …)` floor (`ipsae_v4.py:136`), and that floor — not the clamp — is what makes
`L < 26` return 1.0 (the clamped cubic at `L = 26` is 0.9578). The draft's "within tolerance
everywhere else" also understates the result: the two helpers are **bit-exact** for every
integer `L != 27`, verified by enumeration.

**D5 — Minor: off-by-a-few line ranges.** pDockQ is `:643-670` (draft: `L644-667`); pDockQ2 is
`:672-700` (draft: `L671-698`); LIS is `:702-720` (draft: `L702-718`); ipSAE_d0res starts at
`:786` where `n0res_byres` is computed (draft: `L787`). In each case the draft's range stops
before the zero-contact `else` branch it separately relies on.

**D6 — Not in the draft but worth recording:** the CLAUDE.md pDockQ2 defect is not only the
sigmoid constants. Its `plddt_contacts` construction
(`np.repeat(plddt_A, in_contact.sum(axis=1))` concatenated with
`np.tile(plddt_B, (len(plddt_A),1))[in_contact]`) weights each residue's pLDDT by its number
of contacts, whereas `ipsae.py` takes an unweighted mean over the set of unique interface
residues (`ipsae_v4.py:692`). Fixing only the constants would still leave the wrong `x`.

**D7 — Not in the draft:** CLAUDE.md's `compute_ipsae` is declared to "return a flat dict"
in the surrounding prose but is written as a generator (`yield mode, …`) returning only the
three maxima, with no per-residue arrays, no `d0` values and no `n0` counts. Anything in the
notebook that consumes per-residue ipSAE profiles cannot be built from it as written.

No claim in the draft was found to be *wrong* in the sense of misreporting a formula. The
corrections above are one substantive omission (D1), three citation defects (D2, D3, D5), and
three completeness gaps (D4, D6, D7).
