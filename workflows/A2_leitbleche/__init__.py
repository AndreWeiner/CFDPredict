"""A2_leitbleche as an importable package.

This workflow is normally run as a script (the portal/CLI executes
py_A2_leitbleche.py with the series dir on sys.path), and that keeps working
unchanged. The package form exists so OTHER code can import the concentric
elbow generator without colliding with the identically-named modules that
A2_brand_topology ships -- both trees have gen_blockmesh.py, build_case.py,
params.py and so on, so a flat `import gen_blockmesh` resolves to whichever
happens to be on sys.path first.

Intra-package imports are therefore written as "relative first, absolute as a
fallback": the relative form wins when we are imported as A2_leitbleche.<mod>,
the absolute one when the directory is simply on sys.path (script mode).
"""
