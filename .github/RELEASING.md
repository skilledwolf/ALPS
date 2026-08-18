# Releasing pyalps

1. Update `project.version` in `pyproject.toml` and merge the change through CI.
2. Create a GitHub release whose tag is exactly `v<project.version>`.
3. Publish the GitHub release. Draft releases do not publish packages.
4. Approve the protected `pypi` environment deployment after the source
   distribution and all wheels pass verification.

The release workflow builds one source distribution, then builds every wheel
from that exact archive. Manual dispatch is build-only and can be used to test
release packaging without obtaining a PyPI identity token.

PyPI Trusted Publishing must be restricted to repository `skilledwolf/ALPS`,
workflow `build_wheels.yml`, and environment `pypi`. The publish action emits
PEP 740 attestations by default.
