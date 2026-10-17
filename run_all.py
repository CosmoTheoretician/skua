import subprocess, re, sys
models = sys.argv[1:] or ["lcdm", "olcdm", "wcdm", "cpl", "matter_curv", "signed_x", "dgp", "decay_dm", "spline_x", "janus"]
rows = []
for m in models:
    out = subprocess.run(["python", "test_model.py", m], capture_output=True, text=True).stdout
    print(out)
    c = float(re.search(r"chi2 = ([\d.]+)", out).group(1)); k = int(re.search(r"k = (\d+)", out).group(1))
    fit = out.split("best fit: ")[1].strip()
    rows.append((c + 2 * k, m, c, k, fit))
rows.sort()
lines = ["| model | chi2 | k | AIC | dAIC vs best | best fit |", "|---|---|---|---|---|---|"]
for aic, m, c, k, fit in rows:
    lines.append(f"| {m} | {c:.1f} | {k} | {aic:.1f} | {aic - rows[0][0]:+.1f} | {fit} |")
text = "\n".join(lines); print(text)
open("results/sweep.md", "w").write(text + "\n")
