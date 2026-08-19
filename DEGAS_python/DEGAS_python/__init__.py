from .options import *
from .preprocess import preprocess_counts, normalize_scale


def __getattr__(name):
    """Load torch-dependent public objects only when they are requested.

    This keeps preprocessing and patient-calibration utilities usable on CPU
    analysis nodes where PyTorch is intentionally absent, while retaining the
    historical ``DEGAS_python.run_model`` and ``bagging_all_results`` API.
    """
    if name in {"run_model", "bagging_all_results"}:
        from . import run_module

        return getattr(run_module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
