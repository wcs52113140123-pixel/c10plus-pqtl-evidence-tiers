
import sys; sys.path.insert(0,'.')
from c10common import *
import statsmodels.formula.api as smf
M=derive(load()); D=M[M.v4&~M.v3].dropna(subset=['mlp','absb','maf','logn95']).copy()
for f in ["rep_decode ~ PAV + cr(mlp, df=4) + cr(absb, df=4) + cr(maf, df=4) + cr(logn95, df=4) + decode_has_cis",
          "rep_decode ~ PAV + cr(mlp, df=4) + cr(absb, df=4) + cr(maf, df=4) + cr(logn95, df=4)",
          "rep_decode ~ PAV + cr(mlp, df=4) + cr(np.log(absb), df=4) + cr(maf, df=4) + cr(logn95, df=4) + decode_has_cis"]:
    m=smf.glm(f,D,family=sm.families.Binomial()).fit()
    print(f[:60], m.params['PAV'], m.converged, D.decode_has_cis.value_counts().to_dict(), pd.crosstab(D.decode_has_cis,D.rep_decode).values.tolist())
