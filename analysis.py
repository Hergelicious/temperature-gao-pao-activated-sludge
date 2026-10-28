"""
Temperature and the GAO:PAO balance in full-scale activated sludge.

Reproduces all statistical results reported in the manuscript from the raw
supplementary tables in data/raw/. Writes result tables to results/.

Usage:  python analysis.py
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

MIDAS = RAW / "media-2.xlsx"                            # MiDAS catalog, Supplementary Tables
DK_MAGS = RAW / "41467_2021_22203_MOESM5_ESM.xlsx"     # Singleton et al. 2021, Supplementary Data 3
DK_ABUND = RAW / "41467_2021_22203_MOESM15_ESM.xlsx"   # Singleton et al. 2021, Supplementary Data 13

missing = [f.name for f in (MIDAS, DK_MAGS, DK_ABUND) if not f.exists()]
if missing:
    raise SystemExit(
        "Missing input file(s) in data/raw/: " + ", ".join(missing)
        + "\nDownload them as described in data/raw/README.md and README.md (section 'Input data')."
    )

# ---------------------------------------------------------------------------
# Fixed guild definitions
# ---------------------------------------------------------------------------
PAO = ["Accumulibacter", "Azonexus", "Phosphoribacter", "Lutibacillus"]
GAO = ["Competibacter", "Competibacter_A", "Contendobacter", "Defluviicoccus", "Propionivibrio"]
NOB = ["Nitrospira_A", "Nitrospira_D", "Nitrotoga"]
STORAGE_THRESHOLD = 0.70
X = ["T10", "EBPR", "Nrem", "AS", "lind"]
TEMP_IQR_10C = 1.213  # interquartile temperature range (12.13 degC) in units of 10 degC

# Danish (MiDAS 3.7) equivalents
DK_PAO = ["Ca_Accumulibacter", "Tetrasphaera", "Dechloromonas"]
DK_GAO = ["Ca_Competibacter", "Ca_Contendobacter", "Defluviicoccus", "Propionivibrio"]

# ---------------------------------------------------------------------------
# Global data
# ---------------------------------------------------------------------------
s1 = pd.read_excel(MIDAS, sheet_name="S1")
s1.columns = [c.strip() for c in s1.columns]
s1 = s1[s1.WWTP_id.notna()].set_index("WWTP_id")

s2 = pd.read_excel(MIDAS, sheet_name="S2", header=1)
s2["SA_ID"] = s2["SA_ID"].astype(str).str.strip()
seq_yield = s2.groupby("SA_ID")["Yield (Gb); Len>200bp && >Q10"].sum()

s13 = pd.read_excel(MIDAS, sheet_name="S13")
TRAITS = list(s13.columns[4:])
trait_table = s13.set_index("Genus")[TRAITS] / 100


def genus_table(sheets):
    """Summed relative abundance (%) of MAGs per genus (exact GTDB r226 names) per plant."""
    a = pd.concat(
        [pd.read_excel(MIDAS, sheet_name=s, usecols=["bin", "Relative_abund (%)", "tax_gtdb.r226"]) for s in sheets]
    )
    a["plant"] = a.bin.str.split(".bin").str[0].str.split("-", n=1).str[1]
    a["genus"] = a["tax_gtdb.r226"].str.extract(r"g__([^;]*)")[0]
    return a.groupby(["plant", "genus"])["Relative_abund (%)"].sum().unstack(fill_value=0)


G_ALL = genus_table(["S4", "S5"])   # HQ + MQ MAGs
G_HQ = genus_table(["S4"])          # HQ MAGs only


def storage_genera(thr):
    t = trait_table
    return [g for g in t.index if t.loc[g, "PHA accumulation"] >= thr and t.loc[g, "Glycogen biosynthesis"] >= thr]


def clr_transform(G):
    P = G.div(G.sum(axis=1), axis=0)
    pc = P[P > 0].min().min() / 2
    L = np.log(P + pc)
    return L.sub(L.mean(axis=1), axis=0)


def outcomes(G, thr=STORAGE_THRESHOLD, drop=None):
    """Outcome variables per plant. `drop` removes one genus (leave-one-genus-out)."""
    G = G.drop(columns=[drop], errors="ignore")
    core = G[[g for g in trait_table.index if g in G.columns]]
    core = core.div(core.sum(axis=1), axis=0)
    clr = clr_transform(G)
    o = pd.DataFrame(index=G.index)

    st = [g for g in storage_genera(thr) if g in core.columns]
    s = core[st].sum(axis=1)
    o["H1a"] = np.log((s + 1e-4) / (1 - s + 1e-4))
    o["H1b"] = np.log(G[[g for g in NOB if g in G.columns]].sum(axis=1) + 1e-3)
    ga = [g for g in GAO if g in clr.columns]
    pa = [g for g in PAO if g in clr.columns]
    o["H2"] = clr[ga].mean(axis=1) - clr[pa].mean(axis=1)
    nitrospira = G[[g for g in ["Nitrospira_A", "Nitrospira_D"] if g in G.columns]].sum(axis=1)
    o["H3"] = np.log((G.get("Nitrotoga", 0) + 1e-3) / (nitrospira + 1e-3))

    # alternative closures (sensitivity S2)
    o["H1a_clr"] = clr[[g for g in st if g in clr.columns]].mean(axis=1) - clr[[g for g in clr.columns if g not in st]].mean(axis=1)
    o["H2_core"] = np.log(
        (core[[g for g in ga if g in core.columns]].sum(axis=1) + 1e-4)
        / (core[[g for g in pa if g in core.columns]].sum(axis=1) + 1e-4)
    )
    o["H3_clr"] = clr["Nitrotoga"] - clr[["Nitrospira_A", "Nitrospira_D"]].mean(axis=1)
    o["nitrotoga_present"] = (G.get("Nitrotoga", 0) > 0).astype(int)

    cw = core.dot(trait_table.loc[core.columns])
    cw.columns = ["trait_" + t.replace(" ", "_").replace("/", "_").replace("-", "_") for t in TRAITS]
    return o.join(cw)


def add_covariates(df):
    d = df.join(s1[["Continent", "Country/Region", "Temp_degC", "ClimateZone", "Process_type", "Plant_type", "ind_merged"]])
    d = d[d.Process_type.notna()].copy()
    d["T10"] = d.Temp_degC / 10
    d["EBPR"] = d.Process_type.str.split(",").apply(lambda l: "P" in [x.strip() for x in l]).astype(int)
    d["Nrem"] = d.Process_type.str.contains("N|Anammox").astype(int)
    d["AS"] = (d.Plant_type == "Activated sludge").astype(int)
    d["lind"] = np.log1p(d.ind_merged)
    d["hot_climate"] = d.ClimateZone.astype(str).str[0].isin(["A", "B"]).astype(int)
    d["log_seq_yield"] = np.log(seq_yield.reindex(d.index))
    return d


def fit(y, d, xs=X, term="T10", cluster=False):
    """OLS with HC3 (or country-clustered) standard errors."""
    dd = d.dropna(subset=xs + [y]).copy()
    model = smf.ols(f"{y} ~ " + " + ".join(xs), data=dd)
    if cluster:
        r = model.fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dd["Country/Region"])[0]})
    else:
        r = model.fit(cov_type="HC3")
    ci95 = r.conf_int(0.05).loc[term].values
    ci90 = r.conf_int(0.10).loc[term].values
    return dict(beta=r.params[term], lo95=ci95[0], hi95=ci95[1], lo90=ci90[0], hi90=ci90[1],
                p=r.pvalues[term], se=r.bse[term], n=int(r.nobs))


D = add_covariates(outcomes(G_ALL))

# ---------------------------------------------------------------------------
# Danish temporal baseline (Delta)
# ---------------------------------------------------------------------------
mags = pd.read_excel(DK_MAGS, header=1)
mags["g_midas"] = mags.MiDAS3_7Tax.str.extract(r"g:([^,;]*)")[0]
mags["g_gtdb"] = mags.GTDBTax.str.extract(r"g__([^;]*)")[0]
mags["key"] = mags.MAG.str.replace(r"\.fa$", "", regex=True).str.replace(r"(_[a-z]+)+$", "", regex=True)
ab = pd.read_excel(DK_ABUND, header=1)
ab["key"] = ab.Genome.str.replace(r"(_[a-z]+)+$", "", regex=True)
z = ab.merge(mags[["key", "g_midas", "g_gtdb"]], on="key")
q3 = [c for c in ab.columns if "-Q3-" in str(c)]

dk_m = z.groupby("g_midas")[q3].sum().T
dk_g = z.groupby("g_gtdb")[q3].sum().T
clr_dk = clr_transform(dk_m)
dk = pd.DataFrame(index=q3)
dk["H2"] = clr_dk[[g for g in DK_GAO if g in clr_dk.columns]].mean(axis=1) - clr_dk[[g for g in DK_PAO if g in clr_dk.columns]].mean(axis=1)
dk["H1b"] = np.log(dk_m[["Nitrotoga", "Nitrospira"]].sum(axis=1) + 1e-3)
dk["H3"] = np.log((dk_m["Nitrotoga"] + 1e-3) / (dk_m["Nitrospira"] + 1e-3))
core_dk = dk_g[[g for g in trait_table.index if g in dk_g.columns]]
core_dk = core_dk.div(core_dk.sum(axis=1), axis=0)
st_dk = [g for g in storage_genera(STORAGE_THRESHOLD) if g in core_dk.columns]
s = core_dk[st_dk].sum(axis=1)
dk["H1a"] = np.log((s + 1e-4) / (1 - s + 1e-4))
dk["plant"] = [c.split("_")[0] for c in q3]
dk["year"] = ["20" + c.split("_")[1][:2] for c in q3]
DELTA = dk.groupby("plant")[["H1a", "H1b", "H2", "H3"]].std().median()

dk.rename_axis("sample").reset_index()[["sample", "plant", "year", "H1a", "H1b", "H2", "H3"]].to_csv(
    OUT / "danish_baseline_samples.csv", index=False)
pd.DataFrame({"outcome": DELTA.index, "delta_median_within_plant_SD": DELTA.values,
              "storage_genera_matched": len(st_dk)}).to_csv(OUT / "danish_baseline_delta.csv", index=False)

# ---------------------------------------------------------------------------
# Primary confirmatory models (H1-H3), equivalence, Holm correction
# ---------------------------------------------------------------------------
rows = []
for h in ["H1a", "H1b", "H2", "H3"]:
    r = fit(h, D)
    r["hypothesis"] = h
    if h in ("H1a", "H1b"):
        dl = DELTA[h]
        r["delta"] = dl
        r["equivalent_90CI_within_delta"] = bool(r["lo90"] > -dl and r["hi90"] < dl)
        r["p_TOST"] = max(1 - norm.cdf((r["beta"] + dl) / r["se"]), 1 - norm.cdf((dl - r["beta"]) / r["se"]))
        r["p_confirmatory"] = r["p_TOST"]
    else:
        r["p_confirmatory"] = r["p"]
    rows.append(r)
primary = pd.DataFrame(rows)
primary["p_holm"] = multipletests(primary.p_confirmatory, method="holm")[1]
cols = ["hypothesis", "n", "beta", "lo95", "hi95", "lo90", "hi90", "p", "delta",
        "equivalent_90CI_within_delta", "p_TOST", "p_holm"]
primary[cols].to_csv(OUT / "primary_results.csv", index=False)

# H4: shift across the interquartile temperature range relative to Delta
h4 = []
for h in ["H2", "H3"]:
    r = primary.set_index("hypothesis").loc[h]
    lo, hi = sorted([abs(r.lo95) * TEMP_IQR_10C, abs(r.hi95) * TEMP_IQR_10C])
    h4.append(dict(hypothesis=h, shift_across_IQR=abs(r.beta) * TEMP_IQR_10C, lo95=lo, hi95=hi,
                   delta=DELTA[h], ratio_to_delta=abs(r.beta) * TEMP_IQR_10C / DELTA[h]))
pd.DataFrame(h4).to_csv(OUT / "h4_shift_vs_yearly_variation.csv", index=False)

# H5: EBPR subset, interaction, post hoc power
x_sub = ["T10", "Nrem", "AS", "lind"]
ebpr = fit("H2", D[D.EBPR == 1], x_sub)
no_ebpr = fit("H2", D[D.EBPR == 0], x_sub)
dd = D.dropna(subset=X + ["H2"])
inter = smf.ols("H2 ~ T10*EBPR + Nrem + AS + lind", data=dd).fit(cov_type="HC3")
b_full = primary.set_index("hypothesis").loc["H2", "beta"]
z_a, z_b = norm.ppf(0.975), norm.ppf(0.80)
power = norm.cdf(b_full / ebpr["se"] - z_a) + norm.cdf(-b_full / ebpr["se"] - z_a)
h5 = pd.DataFrame([
    dict(analysis="EBPR plants only", n=ebpr["n"], beta=ebpr["beta"], lo95=ebpr["lo95"], hi95=ebpr["hi95"], p=ebpr["p"]),
    dict(analysis="Plants without EBPR", n=no_ebpr["n"], beta=no_ebpr["beta"], lo95=no_ebpr["lo95"], hi95=no_ebpr["hi95"], p=no_ebpr["p"]),
    dict(analysis="T10 x EBPR interaction (all plants)", n=int(inter.nobs), beta=inter.params["T10:EBPR"],
         lo95=inter.conf_int().loc["T10:EBPR", 0], hi95=inter.conf_int().loc["T10:EBPR", 1], p=inter.pvalues["T10:EBPR"]),
    dict(analysis="Post hoc power of EBPR model to detect full-cohort effect", n=ebpr["n"], beta=power),
    dict(analysis="Minimum detectable effect at 80% power (EBPR model)", n=ebpr["n"], beta=(z_a + z_b) * ebpr["se"]),
])
h5.to_csv(OUT / "h5_ebpr_subset.csv", index=False)

# ---------------------------------------------------------------------------
# Sensitivity analyses S1-S9
# ---------------------------------------------------------------------------
sens, sens_detail = [], []


def record(h, name, r, note=""):
    sens.append(dict(hypothesis=h, analysis=name, n=r.get("n"), beta=r["beta"], lo95=r["lo95"], hi95=r["hi95"], note=note))


for h in ["H1a", "H1b", "H2", "H3"]:
    record(h, "Primary", fit(h, D))
    if h == "H1a":
        for thr in (0.5, 0.9):
            record(h, f"S1 storage threshold {int(thr * 100)}%", fit(h, add_covariates(outcomes(G_ALL, thr))))
    alt = {"H1a": "H1a_clr", "H2": "H2_core", "H3": "H3_clr"}.get(h)
    if alt:
        record(h, "S2 alternative closure", fit(alt, D))
    record(h, "S3 HQ MAGs only", fit(h, add_covariates(outcomes(G_HQ))))
    record(h, "S4 activated sludge plants only", fit(h, D[D.AS == 1], ["T10", "EBPR", "Nrem", "lind"]))
    drops = []
    for c in D.Continent.dropna().unique():
        r = fit(h, D[D.Continent != c])
        drops.append(r)
        sens_detail.append(dict(hypothesis=h, analysis="S5 leave one continent out", left_out=c, **r))
    weakest = min(drops, key=lambda r: abs(r["beta"]))
    record(h, "S5 leave one continent out", weakest,
           f"beta range {min(r['beta'] for r in drops):.2f} to {max(r['beta'] for r in drops):.2f}; row shows weakest")
    if h in ("H1a", "H2"):
        members = [g for g in storage_genera(STORAGE_THRESHOLD) if g in G_ALL.columns] if h == "H1a" \
            else [g for g in PAO + GAO if g in G_ALL.columns]
        drops = []
        for g in members:
            r = fit(h, add_covariates(outcomes(G_ALL, drop=g)))
            drops.append(r)
            sens_detail.append(dict(hypothesis=h, analysis="S6 leave one genus out", left_out=g, **r))
        weakest = min(drops, key=lambda r: abs(r["beta"]))
        record(h, "S6 leave one genus out", weakest,
               f"beta range {min(r['beta'] for r in drops):.2f} to {max(r['beta'] for r in drops):.2f}; row shows weakest")
    record(h, "S7 Koppen A/B vs C/D", fit(h, D, ["hot_climate", "EBPR", "Nrem", "AS", "lind"], term="hot_climate"),
           "contrast, not per 10 degC")
    record(h, "S8 country-clustered SE", fit(h, D, cluster=True))

dd = D.dropna(subset=X + ["log_seq_yield", "nitrotoga_present"])
lg = smf.logit("nitrotoga_present ~ T10 + EBPR + Nrem + AS + lind + log_seq_yield", data=dd).fit(disp=0)
ci = lg.conf_int().loc["T10"].values
record("H3", "S9 presence model with sequencing yield",
       dict(beta=lg.params["T10"], lo95=ci[0], hi95=ci[1], n=int(lg.nobs)),
       f"log-odds; odds ratio per 10 degC {np.exp(lg.params['T10']):.2f} ({np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f})")

pd.DataFrame(sens).to_csv(OUT / "sensitivity_analyses.csv", index=False)
pd.DataFrame(sens_detail).drop(columns=["lo90", "hi90", "se"]).to_csv(OUT / "sensitivity_leave_one_out_detail.csv", index=False)

# ---------------------------------------------------------------------------
# Exploratory trait scan (community-weighted traits, BH-FDR)
# ---------------------------------------------------------------------------
trait_cols = [c for c in D.columns if c.startswith("trait_") and D[c].std() > 1e-6]
rows = []
for t in trait_cols:
    D[t + "_z"] = (D[t] - D[t].mean()) / D[t].std()
    r = fit(t + "_z", D)
    rows.append(dict(trait=t.replace("trait_", "").replace("_", " "), n=r["n"], beta_SD_per_10C=r["beta"],
                     lo95=r["lo95"], hi95=r["hi95"], p=r["p"]))
ts = pd.DataFrame(rows)
ts["q_BH"] = multipletests(ts.p, method="fdr_bh")[1]
ts.sort_values("p").to_csv(OUT / "trait_scan.csv", index=False)

# ---------------------------------------------------------------------------
# Genus-level temperature effects (CLR, all genera)
# ---------------------------------------------------------------------------
clr_all = clr_transform(G_ALL)
core = G_ALL[[g for g in trait_table.index if g in G_ALL.columns]]
core = core.div(core.sum(axis=1), axis=0)
st = [g for g in storage_genera(STORAGE_THRESHOLD) if g in core.columns]
genera = sorted(set([g for g in st if core[g].mean() >= 0.002] + [g for g in PAO + GAO if g in G_ALL.columns])) + NOB
rows = []
for g in genera:
    tmp = D.copy()
    tmp["y"] = clr_all[g].reindex(tmp.index)
    r = fit("y", tmp)
    group = "PAO" if g in PAO else "GAO" if g in GAO else "NOB" if g in NOB else "Other storer"
    rows.append(dict(genus=g, group=group, n=r["n"], beta_CLR_per_10C=r["beta"], lo95=r["lo95"], hi95=r["hi95"],
                     p=r["p"], prevalence=(G_ALL[g] > 0).mean()))
pd.DataFrame(rows).to_csv(OUT / "genus_slopes.csv", index=False)

# ---------------------------------------------------------------------------
# Plant-level analysis dataset
# ---------------------------------------------------------------------------
keep = ["Continent", "Country/Region", "Temp_degC", "ClimateZone", "Process_type", "Plant_type", "ind_merged",
        "EBPR", "Nrem", "AS", "H1a", "H1b", "H2", "H3", "nitrotoga_present", "log_seq_yield"] + trait_cols
plant = D[keep].copy()
plant.insert(0, "complete_case", D[X + ["H2"]].notna().all(axis=1))
plant.rename_axis("WWTP_id").to_csv(OUT / "plant_dataset.csv")

# ---------------------------------------------------------------------------
print("Primary results")
print(primary[["hypothesis", "n", "beta", "lo95", "hi95", "p", "p_holm"]].round(4).to_string(index=False))
print("\nDelta:", DELTA.round(3).to_dict())
print("H5:", h5[["analysis", "beta"]].round(3).to_string(index=False))
print(f"\nResults written to {OUT}")
