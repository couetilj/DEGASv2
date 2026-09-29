from setuptools import setup, find_packages

setup(
    name="DEGAS_python",
    version="0.1",
    packages=find_packages(),
    extras_require={
        "validation": ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "torch", "tqdm"],
        "explain": ["numpy", "torch", "shap"],
    },
    install_requires=[
    ],
)

