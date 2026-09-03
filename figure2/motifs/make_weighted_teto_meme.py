#!/usr/bin/env python3
"""make_weighted_teto_meme.py

Build a weighted PWM MEME motif from a strict MEME motif and a positional weights TSV.

Weighting rule:
  For each position i, new_row = old_row ** w_i; then renormalize to sum=1.
This is equivalent to scaling log-likelihood ratios by w_i when background is uniform.

Inputs:
  --strict-meme   strict motif MEME file (ALPHABET=ACGT; letter-prob matrix)
  --weights-tsv   TSV with columns: pos, base, weight (pos is 1-based)
  --out-meme      output MEME file

Example:
  ./make_weighted_teto_meme.py --strict-meme motifs/tetO_TRE3Gs19.meme \
      --weights-tsv motifs/tetO_positional_weights_19.tsv \
      --out-meme motifs/tetO_TRE3Gs19_weighted.meme
"""
import argparse, re, sys
import pandas as pd
import numpy as np

def parse_matrix(lines, motif_name=None):
    alphabet=None
    for line in lines:
        if line.startswith("ALPHABET="):
            alphabet=line.split("=",1)[1].strip()
    order=list(alphabet) if alphabet else ["A","C","G","T"]

    start=None
    for i,line in enumerate(lines):
        if line.startswith("MOTIF "):
            parts=line.split()
            if motif_name is None or (len(parts)>=2 and parts[1]==motif_name):
                start=i
                break
    if start is None:
        raise SystemExit(f"Motif '{motif_name}' not found")

    mat_header=None
    for j in range(start, len(lines)):
        if lines[j].startswith("letter-probability matrix"):
            mat_header=j
            break
    if mat_header is None:
        raise SystemExit("letter-probability matrix not found")

    m=re.search(r"w=\s*(\d+)", lines[mat_header])
    w=int(m.group(1)) if m else None
    mat=[]
    for k in range(mat_header+1, mat_header+1+w):
        mat.append([float(x) for x in lines[k].strip().split()])
    return order, start, mat_header, np.array(mat), w

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--strict-meme", required=True)
    ap.add_argument("--weights-tsv", required=True)
    ap.add_argument("--out-meme", required=True)
    ap.add_argument("--strict-motif-name", default=None, help="Default: first motif in file")
    ap.add_argument("--out-motif-name", default="tetO_weighted")
    args=ap.parse_args()

    lines=open(args.strict_meme).read().splitlines()
    order, motif_idx, mat_header_idx, mat, w = parse_matrix(lines, args.strict_motif_name)
    weights=pd.read_csv(args.weights_tsv, sep="\t")
    weights=weights.sort_values("pos")
    if len(weights)!=w:
        raise SystemExit(f"weights length {len(weights)} != motif width {w}")
    wvec=weights["weight"].astype(float).to_numpy()

    mat_w=[]
    for i in range(w):
        row=mat[i] ** wvec[i]
        row=row/row.sum()
        mat_w.append(row)
    mat_w=np.array(mat_w)

    out_lines=lines.copy()
    # replace motif line
    out_lines[motif_idx]=f"MOTIF {args.out_motif_name}"
    # replace matrix rows
    for i in range(w):
        out_lines[mat_header_idx+1+i]=" ".join(f"{x:.8f}" for x in mat_w[i])
    open(args.out_meme,"w").write("\n".join(out_lines)+"\n")

if __name__=="__main__":
    main()
