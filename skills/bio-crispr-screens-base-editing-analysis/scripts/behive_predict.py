"""Predict base-editing outcomes for one spacer with BE-Hive (Arbab 2020).

Purpose: BE-Hive needs a fixed 50nt substrate (19nt upstream + 20nt spacer + 3nt PAM + 8nt downstream);
this builds it, checks the spacer offset, then prints total predicted probability and top outcomes.
Inputs:  spacer (20nt), PAM (3nt), 19nt upstream genomic context, 8nt downstream genomic context,
         editor (BE4, ...), cell type (mES, HEK293, U2OS, ...), parent dir of the be_predict_bystander clone.
Usage:   python behive_predict.py --behive-parent /path/to/dl --spacer TGATCACGTAGCATGCACGT --pam TGG \
             --upstream ATGCATGGATCGTAGCTAG --downstream CATGCTAG --editor BE4 --celltype mES
Needs:   the BE-Hive clone (maxwshen/be_predict_bystander) and its Python env (torch, pandas).
"""
import argparse
import sys

ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
ap.add_argument('--behive-parent', required=True, help='parent dir of the cloned be_predict_bystander repo')
ap.add_argument('--spacer', default='TGATCACGTAGCATGCACGT', help='20nt spacer')
ap.add_argument('--pam', default='TGG')
ap.add_argument('--upstream', default='ATGCATGGATCGTAGCTAG', help='19nt of real genomic context immediately 5\' of the spacer')
ap.add_argument('--downstream', default='CATGCTAG', help='8nt of real genomic context immediately 3\' of the PAM')
ap.add_argument('--editor', default='BE4')
ap.add_argument('--celltype', default='mES', help="one of 'mES','HEK293','U2OS',...")
a = ap.parse_args()

sys.path.append(a.behive_parent)  # parent dir of the cloned repo
from be_predict_bystander import predict as bystander_model

spacer = a.spacer  # 20nt
pam = a.pam
upstream_19nt = a.upstream    # 19nt of real genomic context immediately 5' of the spacer
downstream_8nt = a.downstream # 8nt of real genomic context immediately 3' of the PAM
substrate = upstream_19nt + spacer + pam + downstream_8nt
assert len(substrate) == 50

bystander_model.init_model(base_editor=a.editor, celltype=a.celltype)
pred_df, stats = bystander_model.predict(substrate)

# Always cross-check BE-Hive's own read-back against the intended spacer before trusting
# pred_df's position-labeled columns (e.g. 'C4', 'C6') -- a wrong substrate length or
# offset produces a plausible-looking but silently mis-positioned prediction.
assert substrate[19:39] == spacer, "substrate/spacer offset is wrong -- check upstream context length"
print(stats["Total predicted probability"])
print(pred_df.sort_values("Predicted frequency", ascending=False).head(10))
