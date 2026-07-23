import os, re, sys, shutil, tempfile

src = "setupTemplates/otr/defaultCluster/slurm/clusterDict"
text = open(src).read()

repo = "/REPO"; pvbatch = "/PV/bin/pvbatch"; of_root = "/tmp/OF"
solve_prec = "SP"; first_ver = "v2406"
zero_tpl = f"{repo}/zeroTemplates"; post_dir = f"{repo}/postUtilities"
OTR_OF_HOME = "~/openFoam/CODEHOST/OFSource/$VERSION"
TOKEN = r"[^ \"';]*"

def c(t): return t.count("RunFunctions")

print("start", c(text))
text = re.sub(r"caseTemplates='[^']*'", f"caseTemplates='{zero_tpl}'", text); print("after caseTemplates", c(text))
text = re.sub(TOKEN + "postUtilities", post_dir, text); print("after postUtilities", c(text))
text = re.sub(TOKEN + "pvbatch", pvbatch, text); print("after pvbatch", c(text))
print("is_otr", OTR_OF_HOME in text)
text = text.replace(OTR_OF_HOME, f"{of_root}/OpenFOAM-$VERSION"); print("after OFhome", c(text))
text = re.sub(r"VERSION=[A-Za-z]*\d{4}", f"VERSION={first_ver}", text); print("after VERSION", c(text))
solve_sp = "TRUE" if solve_prec.upper() == "SP" else "FALSE"
text = re.sub(r'("precision":"SP=)(?:TRUE|FALSE)(;[^"]*?Running Solve Script)', lambda m: m.group(1)+solve_sp+m.group(2), text); print("after solveSP", c(text))
