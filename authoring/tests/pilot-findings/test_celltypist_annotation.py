import pathlib
import runpy
import sys
import tempfile
import types
import unittest


ROOT = pathlib.Path(__file__).parents[3]
SCRIPT = ROOT / "skills" / "bio-single-cell-cell-annotation" / "examples" / "celltypist_annotation.py"


class FakeMatrix:
    def copy(self):
        return FakeMatrix()


class FakeSeries:
    def __gt__(self, _value):
        return self

    def __invert__(self):
        return self

    def value_counts(self):
        return self

    def to_csv(self, _path):
        return None


class FakeObs(dict):
    def __getitem__(self, key):
        if isinstance(key, FakeSeries):
            return []
        return super().__getitem__(key)


class FakeAdata:
    def __init__(self):
        self.X = FakeMatrix()
        self.layers = {"counts": FakeMatrix()}
        self.obsm = {"X_umap": [[0.0, 0.0], [1.0, 1.0]]}
        self.obs = FakeObs(
            leiden=FakeSeries(),
            majority_voting=FakeSeries(),
            conf_score=FakeSeries(),
        )

    def __len__(self):
        return 2

    def write(self, _path):
        return None


class Predictions:
    def __init__(self, adata):
        self.adata = adata

    def to_adata(self):
        return self.adata


class CellTypistContractTests(unittest.TestCase):
    def run_example(self, adata):
        calls = {"normalize_total": [], "log1p": [], "annotate": [], "model_load": []}

        scanpy = types.ModuleType("scanpy")
        scanpy.read_h5ad = lambda _path: adata
        scanpy.pp = types.SimpleNamespace(
            normalize_total=lambda value, **kwargs: calls["normalize_total"].append((value, kwargs)),
            log1p=lambda value: calls["log1p"].append(value),
        )
        scanpy.pl = types.SimpleNamespace(umap=lambda *_args, **_kwargs: None)

        def reject_download(**_kwargs):
            self.fail("The example attempted a network model download")

        def annotate(value, **kwargs):
            calls["annotate"].append((value, kwargs))
            return Predictions(value)

        celltypist = types.ModuleType("celltypist")

        pyplot = types.ModuleType("matplotlib.pyplot")
        pyplot.subplots = lambda *_args, **_kwargs: (object(), [object(), object()])
        pyplot.tight_layout = lambda: None
        pyplot.savefig = lambda *_args, **_kwargs: None
        pyplot.close = lambda: None
        matplotlib = types.ModuleType("matplotlib")
        matplotlib.pyplot = pyplot

        with tempfile.TemporaryDirectory() as model_cache:
            model_file = pathlib.Path(model_cache) / "Immune_All_Low.pkl"
            model_file.touch()
            celltypist.models = types.SimpleNamespace(
                models_path=model_cache,
                download_models=reject_download,
                Model=types.SimpleNamespace(
                    load=lambda **kwargs: calls["model_load"].append(kwargs) or object()
                ),
            )
            celltypist.annotate = annotate

            replacements = {
                "scanpy": scanpy,
                "celltypist": celltypist,
                "matplotlib": matplotlib,
                "matplotlib.pyplot": pyplot,
            }
            previous = {name: sys.modules.get(name) for name in replacements}
            old_argv = sys.argv
            try:
                sys.modules.update(replacements)
                sys.argv = [str(SCRIPT)]
                runpy.run_path(str(SCRIPT), run_name="__main__")
            finally:
                sys.argv = old_argv
                for name, value in previous.items():
                    if value is None:
                        sys.modules.pop(name, None)
                    else:
                        sys.modules[name] = value

        return calls

    def test_normalizes_counts_and_uses_seeded_over_clustering_without_download(self):
        adata = FakeAdata()
        calls = self.run_example(adata)

        self.assertEqual(len(calls["normalize_total"]), 1)
        self.assertEqual(calls["normalize_total"][0][1], {"target_sum": 1e4})
        self.assertEqual(calls["log1p"], [adata])
        self.assertEqual(len(calls["annotate"]), 1)
        self.assertEqual(calls["annotate"][0][1]["over_clustering"], "leiden")
        self.assertEqual(len(calls["model_load"]), 1)
        self.assertIn("/", calls["model_load"][0]["model"])

    def test_rejects_missing_required_analysis_state(self):
        cases = {
            "counts layer": (lambda adata: adata.layers.pop("counts"), "raw counts"),
            "seeded clusters": (lambda adata: adata.obs.pop("leiden"), "seeded over-clustering"),
            "UMAP": (lambda adata: adata.obsm.pop("X_umap"), "X_umap"),
        }

        for name, (remove, expected_message) in cases.items():
            with self.subTest(name=name):
                adata = FakeAdata()
                remove(adata)
                with self.assertRaisesRegex(ValueError, expected_message):
                    self.run_example(adata)


if __name__ == "__main__":
    unittest.main()
