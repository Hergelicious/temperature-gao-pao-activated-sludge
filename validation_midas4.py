"""
Independent validation of the temperature effect on the GAO:PAO balance in the
MiDAS 4 global 16S rRNA gene amplicon survey (Dueholm et al. 2022).

Inputs (see data/raw/README.md):
  data/raw/midas4/  extracted contents of Dueholm2021a_data.zip (figshare 16566408)
  data/raw/media-2.xlsx  MiDAS global genome catalog Supplementary Tables (plant IDs of the discovery cohort)
Outputs: results/validation_*.csv
"""
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
M4 = RAW / "midas4"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


def find(pattern):
    hits = sorted(M4.rglob(pattern))
    if not hits:
        raise SystemExit(f"Missing MiDAS 4 file matching '{pattern}' in {M4}. See data/raw/README.md.")
    return hits[0]


SEQMETA = find("*Sequencing_metadata.txt")
PLANTMETA = find("DataS1*.txt")
CATALOG = RAW / "media-2.xlsx"

# MiDAS 4.8 genus names for the fixed guilds (same mapping as the Danish baseline)
PAO = ["Ca_Accumulibacter", "Tetrasphaera", "Dechloromonas"]
GAO = ["Ca_Competibacter", "Ca_Contendobacter", "Defluviicoccus", "Propionivibrio"]
X = ["T10", "EBPR", "Nrem", "AS", "lind"]
MIN_READS = 10000
DISCOVERY_BETA = 1.0456  # full-cohort H2 effect per 10 degC in the genome catalog


# ---------------------------------------------------------------------------
# Load amplicon data and aggregate to genus
# ---------------------------------------------------------------------------
def read_sintax(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            asv = f[0].split(";")[0]
            tax = f[3] if len(f) > 3 and f[3] else ""
            g = re.search(r"g:([^,]+)", tax)
            rows.append((asv, g.group(1) if g else np.nan))
    return pd.Series(dict(rows), name="genus")


def plant_profiles(region):
    """Genus proportions per plant (mean of replicate libraries) for amplicon region 'V13' or 'V4'."""
    tab = pd.read_csv(find(f"{region}_ASVtab.txt"), sep="\t", index_col=0).select_dtypes("number")
    genus = read_sintax(find(f"{region}ASV_vs_MiDAS_4.8.sintax")).reindex(tab.index)
    counts = tab.groupby(genus.fillna("unclassified")).sum().T          # libraries x genera
    counts = counts[counts.sum(axis=1) >= MIN_READS]
    counts = counts.drop(columns=["unclassified"], errors="ignore")
    seq = pd.read_csv(SEQMETA, sep="\t")
    seq = seq[[f"{region}_seq_id", "WWTP_id"]].dropna().set_index(f"{region}_seq_id")["WWTP_id"]
    counts = counts[counts.index.isin(seq.index)]
    P = counts.div(counts.sum(axis=1), axis=0)
    return P.groupby(seq.reindex(P.index)).mean()


P = plant_profiles("V13")


def clr_transform(P):
    P = P.div(P.sum(axis=1), axis=0)
    pc = P[P > 0].min().min() / 2
    L = np.log(P + pc)
    return L.sub(L.mean(axis=1), axis=0)


def balance(clr, gao, pao):
    return clr[[g for g in gao if g in clr]].mean(axis=1) - clr[[g for g in pao if g in clr]].mean(axis=1)


clr = clr_transform(P)
o = pd.DataFrame(index=P.index)
o["H2"] = balance(clr, GAO, PAO)
o["H2_noDechloromonas"] = balance(clr, GAO, ["Ca_Accumulibacter", "Tetrasphaera"])
o["H2_withMicropruina"] = balance(clr, GAO + ["Micropruina"], PAO)
o["H2_Accumulibacter"] = clr["Ca_Competibacter"] - clr["Ca_Accumulibacter"]
o["H3"] = np.log((P.get("Nitrotoga", 0) + 1e-5) / (P.get("Nitrospira", 0) + 1e-5))
for g in PAO + GAO + ["Micropruina", "Nitrotoga", "Nitrospira"]:
    if g in clr:
        o["clr_" + g] = clr[g]

# ---------------------------------------------------------------------------
# Plant metadata and covariates
# ---------------------------------------------------------------------------
meta = pd.read_csv(PLANTMETA, sep="\t")
meta.columns = [c.strip() for c in meta.columns]
meta = meta.drop_duplicates("WWTP_id").set_index("WWTP_id")
d = o.join(meta, how="inner")
d = d[d.Process_type.notna() & d.Temp_degC.notna()].copy()
d["Temp_degC"] = pd.to_numeric(d.Temp_degC, errors="coerce")
d = d[d.Temp_degC.notna()]
d["T10"] = d.Temp_degC / 10
d["EBPR"] = d.Process_type.astype(str).str.split(",").apply(lambda l: "P" in [x.strip() for x in l]).astype(int)
d["Nrem"] = d.Process_type.astype(str).str.contains("N|Anammox").astype(int)
d["AS"] = (d.Plant_type == "Activated sludge").astype(int)
d["lind"] = np.log1p(pd.to_numeric(d.ind_merged.replace({"20-50": "35"}), errors="coerce"))  # range coded as midpoint

cat = pd.read_excel(CATALOG, sheet_name="S1")
cat_ids = set(cat.WWTP_id.dropna())
d["cohort"] = np.where(d.index.isin(cat_ids), "catalog", "independent")
country_col = [c for c in d.columns if c.startswith("Country")][0]


def fit(y, df, xs=X, term="T10", cluster=False, cov="HC3"):
    dd = df.dropna(subset=[x for x in xs if x in df.columns] + [y]).copy()
    m = smf.ols(f"{y} ~ " + " + ".join(xs), data=dd)
    r = m.fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(dd[country_col])[0]}) if cluster else m.fit(cov_type=cov)
    ci = r.conf_int(0.05).loc[term].values
    return dict(beta=r.params[term], lo95=ci[0], hi95=ci[1], p=r.pvalues[term], se=r.bse[term], n=int(r.nobs))


ind = d[d.cohort == "independent"]
rows = []


def rec(name, r, **kw):
    rows.append(dict(analysis=name, **r, **kw))


rec("V1 independent plants, GAO:PAO balance (primary replication)", fit("H2", ind))
rec("V2 independent EBPR plants", fit("H2", ind[ind.EBPR == 1], ["T10", "Nrem", "AS", "lind"]))
rec("V2 independent plants without EBPR", fit("H2", ind[ind.EBPR == 0], ["T10", "Nrem", "AS", "lind"]))
inter = smf.ols("H2 ~ T10*EBPR + Nrem + AS + lind", data=ind.dropna(subset=X + ["H2"])).fit(cov_type="HC3")
rec("V2 T10 x EBPR interaction (independent plants)",
    dict(beta=inter.params["T10:EBPR"], lo95=inter.conf_int().loc["T10:EBPR", 0],
         hi95=inter.conf_int().loc["T10:EBPR", 1], p=inter.pvalues["T10:EBPR"], se=inter.bse["T10:EBPR"], n=int(inter.nobs)))
rec("V2 all amplicon EBPR plants (independent + catalog)", fit("H2", d[d.EBPR == 1], ["T10", "Nrem", "AS", "lind"]))
rec("V3 sensitivity: Dechloromonas excluded from PAOs", fit("H2_noDechloromonas", ind))
rec("V3 sensitivity: Micropruina added to GAOs", fit("H2_withMicropruina", ind))
rec("V3 sensitivity: Competibacter vs Accumulibacter only", fit("H2_Accumulibacter", ind))
rec("V3 sensitivity: activated sludge plants only", fit("H2", ind[ind.AS == 1], ["T10", "EBPR", "Nrem", "lind"]))
rec("V3 sensitivity: country-clustered errors", fit("H2", ind, cluster=True))
rec("V3 sensitivity: no industrial-load covariate", fit("H2", ind, ["T10", "EBPR", "Nrem", "AS"]))
for cont in sorted(ind.Continent.dropna().unique()):
    rec(f"V3 sensitivity: without {cont}", fit("H2", ind[ind.Continent != cont]))
rec("V3 sensitivity: country fixed effects (within-country temperature effect)",
    fit("H2", ind.assign(country=ind[country_col]).groupby("country").filter(lambda g: len(g) >= 3),
        X + ["C(country)"], cov="HC1"))
rec("V3 sensitivity: independent EBPR plants, country-clustered errors",
    fit("H2", ind[ind.EBPR == 1], ["T10", "Nrem", "AS", "lind"], cluster=True))
ind_k = ind.assign(hot_climate=ind.ClimateZone.astype(str).str[0].isin(["A", "B"]).astype(int))
ind_k = ind_k[ind.ClimateZone.notna()]
rec("V3 sensitivity: Koppen A/B vs C/D in place of temperature",
    fit("H2", ind_k, ["hot_climate", "EBPR", "Nrem", "AS", "lind"], term="hot_climate"))
P4 = plant_profiles("V4")
v4 = d.drop(columns=["H2"]).join(balance(clr_transform(P4), GAO, PAO).rename("H2"), how="inner")
rec("V3 sensitivity: V4 primer set instead of V1-V3", fit("H2", v4[v4.cohort == "independent"]))
rec("V4 independent plants, Nitrotoga:Nitrospira log-ratio", fit("H3", ind))
rec("V5 catalog plants re-measured by amplicon sequencing", fit("H2", d[d.cohort == "catalog"]))

# power of the independent EBPR subset to detect the discovery effect
e = fit("H2", ind[ind.EBPR == 1], ["T10", "Nrem", "AS", "lind"])
za, zb = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
power = stats.norm.cdf(DISCOVERY_BETA / e["se"] - za) + stats.norm.cdf(-DISCOVERY_BETA / e["se"] - za)
rows.append(dict(analysis="Power of independent EBPR model to detect discovery effect", beta=power, n=e["n"]))
rows.append(dict(analysis="Minimum detectable effect at 80% power (independent EBPR model)", beta=(za + zb) * e["se"], n=e["n"]))
res = pd.DataFrame(rows)
res.to_csv(OUT / "validation_results.csv", index=False)

# genus-level slopes in independent plants
gs = []
for g in [c for c in d.columns if c.startswith("clr_")]:
    r = fit(g, ind)
    gs.append(dict(genus=g[4:], **r))
pd.DataFrame(gs).to_csv(OUT / "validation_genus_slopes.csv", index=False)

# cross-method concordance in plants shared with the genome catalog
mag = pd.read_csv(OUT / "plant_dataset.csv", index_col=0) if (OUT / "plant_dataset.csv").exists() else None
if mag is not None:
    shared = d[d.cohort == "catalog"][["H2"]].join(mag[["H2"]], rsuffix="_mag", how="inner").dropna()
    rho, prho = stats.spearmanr(shared.H2, shared.H2_mag)
    r, pr = stats.pearsonr(shared.H2, shared.H2_mag)
    pd.DataFrame([dict(n=len(shared), spearman_rho=rho, spearman_p=prho, pearson_r=r, pearson_p=pr)]).to_csv(
        OUT / "validation_method_concordance.csv", index=False)

d.drop(columns=[c for c in meta.columns if c not in
                ["Continent", country_col, "Temp_degC", "Process_type", "Plant_type", "ind_merged", "ClimateZone"]],
       errors="ignore").to_csv(OUT / "validation_plant_dataset.csv")
print(res.to_string())
