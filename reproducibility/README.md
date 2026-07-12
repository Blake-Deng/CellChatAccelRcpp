# Reproducible Minimal Example

This directory provides a small offline example for reviewer checks:

- `CellChatAccelRcpp_reproducible_minimal_example_a2af664.tar.gz`

The archive is pinned to GitHub commit `a2af6647db29bf6d1b59a4f42ead6e848ed9ab8e`
and includes a prepared PBMC3k CellChat object, an offline package source
tarball, original CellChat and `sparse_stream` run scripts, numerical agreement
checks, runtime and peak-RSS records, expected outputs and SHA-256 checksums.

Run from a clean directory:

```bash
tar -xzf CellChatAccelRcpp_reproducible_minimal_example_a2af664.tar.gz
cd CellChatAccelRcpp_reproducible_minimal_example
bash run_example.sh
```

Expected final line:

```text
VALIDATION PASSED
```

Archive SHA-256:

```text
1a6fca64fb83cd138fbdeb2f22d3a7e3f8152439c947581cdbb745bd8a10fe4e  CellChatAccelRcpp_reproducible_minimal_example_a2af664.tar.gz
```

This is a functional installation and numerical-equivalence example, not a
large-scale performance benchmark. Runtime and memory values depend on the
local machine and system load.
