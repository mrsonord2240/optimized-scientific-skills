import pandas as pd
import numpy as np

def load_iqtree_state(state_file):
    '''Validate an IQ-TREE .state table and group rows by node.

    The provider format has Node, Site, State, and one or more p_<state>
    posterior columns. Invalid or incomplete tables fail with a stable error.
    '''
    df = pd.read_csv(state_file, sep='\t', comment='#')
    df.columns = [c.strip() for c in df.columns]
    required = {'Node', 'Site', 'State'}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"IQ-TREE state table missing required columns: {', '.join(missing)}")
    state_cols = [c for c in df.columns if c.startswith('p_')]
    if not state_cols:
        raise ValueError('IQ-TREE state table requires at least one p_<state> posterior column')
    if df.empty:
        raise ValueError('IQ-TREE state table contains no reconstruction rows')
    if df[['Node', 'Site', 'State']].isna().any().any():
        raise ValueError('IQ-TREE state table has missing Node, Site, or State values')
    if (df['Node'].astype(str).str.strip().eq('').any()
            or df['State'].astype(str).str.strip().eq('').any()):
        raise ValueError('IQ-TREE state table has blank Node or State values')
    try:
        df['Site'] = pd.to_numeric(df['Site'], errors='raise')
        df[state_cols] = df[state_cols].apply(pd.to_numeric, errors='raise')
    except (TypeError, ValueError) as exc:
        raise ValueError(f'IQ-TREE state table has non-numeric Site or posterior values: {exc}') from exc
    if not np.isfinite(df[['Site', *state_cols]].to_numpy(dtype=float)).all():
        raise ValueError('IQ-TREE state table contains non-finite Site or posterior values')
    if not np.equal(df['Site'].to_numpy(), np.floor(df['Site'].to_numpy())).all():
        raise ValueError('IQ-TREE Site values must be integers')
    posterior = df[state_cols].to_numpy(dtype=float)
    if ((posterior < 0.0) | (posterior > 1.0)).any():
        raise ValueError('IQ-TREE posterior probabilities must be between 0 and 1')
    row_sums = posterior.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-3, rtol=1e-3):
        bad = int(np.flatnonzero(~np.isclose(row_sums, 1.0, atol=1e-3, rtol=1e-3))[0]) + 1
        raise ValueError(f'IQ-TREE posterior probabilities do not sum to 1 at data row {bad}')
    df['max_post'] = df[state_cols].max(axis=1)
    maxima = df[state_cols].max(axis=1).to_numpy(dtype=float)
    state_values = df['State'].astype(str).to_numpy()
    if any(not np.isclose(df.iloc[i][f'p_{state_values[i]}'], maxima[i], atol=1e-12)
           if f'p_{state_values[i]}' in state_cols else True
           for i in range(len(df))):
        raise ValueError('IQ-TREE State must be a maximum-posterior state in each row')
    return df.groupby('Node')
